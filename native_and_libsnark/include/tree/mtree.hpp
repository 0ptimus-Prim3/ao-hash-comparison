#pragma once

#include "util/const_math.hpp"
#include "util/string_utils.hpp"

#include <array>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>

#if __cplusplus >= 202002L
    #include <ranges>
#endif

template<typename Hash>
class MTreeNode
{
public:
    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;

private:
    std::array<uint8_t, Hash::DIGEST_SIZE> digest;
    MTreeNode *f;
    std::array<MTreeNode *, ARITY> c;
    size_t depth;

    template<size_t, typename>
    friend class MTree;

    template<size_t, typename>
    friend class MTreePath;

public:
    MTreeNode() = default;
    MTreeNode(const MTreeNode &) = default;
    MTreeNode(MTreeNode &&) = default;
    MTreeNode &operator=(const MTreeNode &) = default;
    MTreeNode &operator=(MTreeNode &&) = default;

    MTreeNode(const void *data, size_t depth) : digest{}, f{}, c{}, depth{depth}
    {
        Hash::hash_oneblock(this->digest.data(), data);
    }

    MTreeNode(const std::array<const void *, ARITY> &data, size_t depth) :
        digest{}, f{}, c{}, depth{depth}
    {
        std::array<uint8_t, Hash::BLOCK_SIZE> block;

        for (size_t i = 0; i < ARITY; ++i)
            memcpy(block.data() + i * Hash::DIGEST_SIZE, data[i], Hash::DIGEST_SIZE);

        Hash::hash_oneblock(this->digest.data(), block.data());
    }

    const auto &get_digest() const { return digest; }
    const MTreeNode *parent() const { return f; }
    const MTreeNode *child(size_t i) const { return c[i]; }

    friend std::ostream &operator<<(std::ostream &os, const MTreeNode &node)
    {
        for (size_t i = 0; i < node.depth; ++i)
            os << "    ";

        os << "*: " << hexdump(node.digest) << '\n';

        for (size_t i = 0; i < ARITY; ++i)
            if (node.c[i] != nullptr)
                os << *node.c[i];

        return os;
    }
};

template<size_t height, typename Hash>
class MTree
{
public:
    using Node = MTreeNode<Hash>;

    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;
    static constexpr size_t LEAVES_N = pow(ARITY, height - 1);
    static constexpr size_t NODES_N = pow_sum(ARITY, (size_t)0, height);
    static constexpr size_t INPUT_SIZE = LEAVES_N * Hash::BLOCK_SIZE;

private:
    /*
    Nodes layout is as follows:
    - The first LEAVES_N nodes contain the leaves
    - The remaining nodes are the internal nodes of the tree
    */
    std::vector<Node> nodes;
    Node *root;

public:
    MTree() = default;
    MTree(const MTree &other) : nodes{other.nodes}, root{&nodes.back()}
    {
        // fixup pointers
        ptrdiff_t off = nodes.data() - other.nodes.data();

        for (size_t i = 0; i < NODES_N; ++i)
        {
            if (nodes[i].f)
                nodes[i].f += off;

            for (size_t j = 0; j < ARITY; ++j)
                if (nodes[i].c[j])
                    nodes[i].c[j] += off;
        }
    }
    MTree(MTree &&other) = default;

    MTree &operator=(const MTree &other)
    {
        nodes = other.nodes;
        root = &nodes.back();

        // fixup pointers
        ptrdiff_t off = nodes.data() - other.nodes.data();

        for (size_t i = 0; i < NODES_N; ++i)
        {
            if (nodes[i].f)
                nodes[i].f += off;

            for (size_t j = 0; j < ARITY; ++j)
                if (nodes[i].c[j])
                    nodes[i].c[j] += off;
        }

        return *this;
    }

    MTree &operator=(MTree &&other) = default;


#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    MTree(const Range &range) :
        MTree{std::ranges::cdata(range),
              std::ranges::size(range) * sizeof(*std::ranges::cdata(range))}
    {}
#endif

    template<typename Iter>
    MTree(const Iter begin, const Iter end) :
        MTree{&*begin, std::distance(begin, end) * sizeof(*begin)}
    {}

    MTree(const void *vdata, size_t sz) : nodes(NODES_N), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "MTree: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

#ifdef MULTICORE
    #pragma omp parallel for
#endif
        // add leaves
        for (size_t i = 0; i < LEAVES_N; ++i)
            this->nodes[i] = Node{data + i * Hash::BLOCK_SIZE, depth};

        // build tree bottom-up
        for (size_t i = 0, len = LEAVES_N; depth > 0; i += len * ARITY)
        {
            size_t last = i + len;
            --depth;
            len /= ARITY;

#ifdef MULTICORE
    #pragma omp parallel for
#endif
            for (size_t j = 0; j < len; ++j)
            {
                std::array<const void *, ARITY> children;

                for (size_t k = 0; k < ARITY; ++k)
                    children[k] = this->nodes[i + j * ARITY + k].digest.data();

                this->nodes[last + j] = Node{children, depth};

                for (size_t k = 0; k < ARITY; ++k)
                {
                    this->nodes[last + j].c[k] = &this->nodes[i + j * ARITY + k];
                    this->nodes[i + j * ARITY + k].f = &this->nodes[last + j];
                }
            }
        }
    }

    const uint8_t *digest() const { return root->digest.data(); }

    const Node *get_node(size_t i) const { return &nodes[i]; }

    friend std::ostream &operator<<(std::ostream &os, const MTree &tree)
    {
        if (!tree.root)
            return os << "*:";

        return os << *tree.root;
    }
};


template<size_t height, typename Hash>
class MTreePath
{
public:
    using Node = MTreeNode<Hash>;

    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;
    static constexpr size_t LEAVES_N = pow(ARITY, height - 1);
    static constexpr size_t NODES_N = height;
    static constexpr size_t INPUT_N = (ARITY - 1) * (height - 1);
    static constexpr size_t INPUT_SIZE = INPUT_N * Hash::DIGEST_SIZE + Hash::BLOCK_SIZE;

private:
    /*
    Nodes layout is as follows:
    - nodes contain the path to the root
    */
    std::vector<Node> nodes;
    Node *root;

public:
    MTreePath() = default;

    MTreePath(const MTreePath &other) : nodes{other.nodes}, root{&nodes.back()}
    {
        // fixup pointers
        ptrdiff_t off = nodes.data() - other.nodes.data();

        for (size_t i = 0; i < NODES_N; ++i)
        {
            if (nodes[i].f)
                nodes[i].f += off;

            for (size_t j = 0; j < ARITY; ++j)
                if (nodes[i].c[j])
                    nodes[i].c[j] += off;
        }
    }
    MTreePath(MTreePath &&other) = default;

    MTreePath &operator=(const MTreePath &other)
    {
        nodes = other.nodes;
        root = &nodes.back();

        // fixup pointers
        ptrdiff_t off = nodes.data() - other.nodes.data();

        for (size_t i = 0; i < NODES_N; ++i)
        {
            if (nodes[i].f)
                nodes[i].f += off;

            for (size_t j = 0; j < ARITY; ++j)
                if (nodes[i].c[j])
                    nodes[i].c[j] += off;
        }

        return *this;
    }
    MTreePath &operator=(MTreePath &&other) = default;


#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    MTreePath(const Range &range, size_t idx = 0) :
        MTreePath{std::ranges::cdata(range),
                  std::ranges::size(range) * sizeof(*std::ranges::cdata(range))}
    {}
#endif

    template<typename Iter>
    MTreePath(const Iter begin, const Iter end, size_t idx = 0) :
        MTreePath{&*begin, std::distance(begin, end) * sizeof(*begin)}
    {}

    MTreePath(const void *vdata, size_t sz, size_t idx = 0) : nodes(NODES_N), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "MTreePath: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

        // bootstrap first node of the path
        this->nodes[0] = Node{data, depth};

        // build tree bottom-up
        for (size_t i = 1; i < height; ++i, idx /= ARITY)
        {
            std::array<const void *, ARITY> children;
            size_t off = Hash::BLOCK_SIZE + (i - 1) * (ARITY - 1) * Hash::DIGEST_SIZE;
            size_t j = idx % ARITY;

            --depth;

            // build correct permutation depending on idx
            for (size_t k = 0; k < j; ++k)
                children[k] = data + off + k * Hash::DIGEST_SIZE;

            children[j] = this->nodes[i - 1].digest.data();

            for (size_t k = j + 1; k < ARITY; ++k)
                children[k] = data + off + (k - 1) * Hash::DIGEST_SIZE;

            // create node and link children
            this->nodes[i] = Node{children, depth};

            this->nodes[i].c[j] = &this->nodes[i - 1];
            this->nodes[i - 1].f = &this->nodes[i];
        }
    }

    const uint8_t *digest() const { return root->digest.data(); }

    const Node *get_node(size_t i) const { return &nodes[i]; }

    friend std::ostream &operator<<(std::ostream &os, const MTreePath &tree)
    {
        if (!tree.root)
            return os;

        return os << *tree.root;
    }
};
