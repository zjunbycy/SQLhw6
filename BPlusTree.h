#ifndef BPLUSTREE_H_INCLUDED
#define BPLUSTREE_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <cstdint>
#include <unordered_set>

namespace hw6 {

	// Z-Curve 空间填充曲线工具类
	class SpaceFilling {
	public:
		// 计算 Z-order (Morton) 值
		// 将2D空间坐标(x, y)映射到1D的Z值
		static uint64_t xyToZOrder(uint32_t x, uint32_t y);

		// 将浮点坐标归一化为整数网格坐标
		static uint32_t normalizeCoord(double coord, double minCoord, double maxCoord, uint32_t gridSize);

		// 计算特征的Z值（使用其包围盒中心点）
		static uint64_t computeZValue(const Feature& feature, const Envelope& bbox, uint32_t gridSize);
	};

	// B+树节点
	class BPlusNode {
	private:
		bool isLeaf;
		BPlusNode* parent;
		std::vector<uint64_t> keys;           // Z值作为键
		std::vector<Feature> features;        // 叶节点存储的特征
		std::vector<BPlusNode*> children;     // 内部节点的子节点指针
		BPlusNode* next;                      // 叶节点链表指针

	public:
		BPlusNode(bool leaf = true);
		~BPlusNode();

		bool getIsLeaf() const { return isLeaf; }
		BPlusNode* getParent() const { return parent; }
		void setParent(BPlusNode* p) { parent = p; }
		BPlusNode* getNext() const { return next; }
		void setNext(BPlusNode* n) { next = n; }

		size_t getKeyNum() const { return keys.size(); }
		uint64_t getKey(size_t i) const { return keys[i]; }
		const std::vector<uint64_t>& getKeys() const { return keys; }

		size_t getFeatureNum() const { return features.size(); }
		const Feature& getFeature(size_t i) const { return features[i]; }
		const std::vector<Feature>& getFeatures() const { return features; }

		size_t getChildNum() const { return children.size(); }
		BPlusNode* getChild(size_t i) const { return children[i]; }

		// 新增：用于构建新根节点的公共方法
		void addKey(uint64_t key) { keys.push_back(key); }
		void addChild(BPlusNode* child) { 
			children.push_back(child); 
			if (child) child->setParent(this);
		}

		void insertLeaf(uint64_t key, const Feature& feature, size_t capacity, BPlusNode*& newNode);
		void insertInternal(uint64_t key, BPlusNode* child);
		void split(size_t capacity, BPlusNode*& newNode, uint64_t& midKey);

		// 查找应该插入的子节点
		size_t findChildIndex(uint64_t key) const;
	};

	// B+树空间索引
	class BPlusTree : public Tree {
	private:
		BPlusNode* root;
		uint32_t gridSize;  // 空间网格划分粒度（用于Z-order计算）
		
		// 插入辅助函数
		void insertHelper(uint64_t zValue, const Feature& feature, BPlusNode* node, BPlusNode*& newNode, uint64_t& midKey);

	public:
		BPlusTree();
		BPlusTree(size_t cap, uint32_t grid = 1024);
		~BPlusTree();

		// 设置网格大小
		void setGridSize(uint32_t grid) { gridSize = grid; }

		// 构建B+树索引（基于Z-order排序）
		virtual bool constructTree(const std::vector<Feature>& features) override;

		// 区域查询
		virtual void rangeQuery(const Envelope& rect, std::vector<Feature>& features) override;

		// 最近邻查询
		virtual bool NNQuery(double x, double y, std::vector<Feature>& features) override;

		// 空间连接
		void spatialJoin(const std::vector<Feature>& A,
			const std::vector<Feature>& B,
			double dist,
			std::vector<std::pair<Feature, Feature>>& result);

		// 统计节点
		virtual void countNode(int& interiorNum, int& leafNum) override;
		void countNodeHelper(BPlusNode* node, int& interiorNum, int& leafNum);

		// 计算高度
		virtual void countHeight(int& height) override;
		int countHeightHelper(BPlusNode* node);

		// 绘制（可选）
		virtual void draw() override;

	public:
		static void test(int t);
		static void analyse();
	};

} // namespace hw6

#endif // BPLUSTREE_H_INCLUDED