#ifndef RTREE_H_INCLUDED
#define RTREE_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"

#include <algorithm>
#include <array>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "CMakeIn.h"

namespace hw6 {

	/// <summary>
	/// </summary>
	class RNode {
	private:
		RNode* parent = nullptr;
		int maxChildren;
		std::vector<RNode*> children;
		int childrenNum = 0;
		Envelope bbox;
		std::vector<Feature> features;

	public:
		RNode() = delete;
		RNode(const Envelope& box) : bbox(box) {}

		bool isLeafNode() { return childrenNum == 0; }

		const Envelope& getEnvelope() { return bbox; }

		RNode* getParent() { return parent; }

		void setEnvelope(const Envelope& box) { bbox = box; }

		RNode* getChildNode(size_t i) {
			return i < childrenNum ? children[i] : nullptr;
		}

		const RNode* getChildNode(size_t i) const {
			return i < childrenNum ? children[i] : nullptr;
		}

		int getChildNum() const { return childrenNum; }

		size_t getFeatureNum() const { return features.size(); }

		const Feature& getFeature(size_t i) const { return features[i]; }

		const std::vector<Feature>& getFeatures() const { return features; }

		void add(const Feature& f) { features.push_back(f); }

		void add(RNode* child);

		void remove(const Feature& f);

		void remove(RNode* child);

		Feature popBackFeature();

		RNode* popBackChildNode();

		void countNode(int& interiorNum, int& leafNum);

		int countHeight(int height);

		void draw();

		void rangeQuery(const Envelope& rect, std::vector<Feature>& features);

		RNode* pointInLeafNode(double x, double y);
	};

	class RTree : public Tree {
	private:
		int maxChildren;
		RNode* root = nullptr;

	public:
		// 默认构造函数，maxChildren 默认为 8
		RTree() : RTree(8) {}
		explicit RTree(int maxChildren);
		~RTree() {
			if (root != nullptr) delete root;
			root = nullptr;
		}

		void setCapacity(int capacity) override {
			// DO NOTHING, since capacity is immutable in R tree
		}

		bool constructTree(const std::vector<Feature>& features) override;

		void countNode(int& interiorNum, int& leafNum) override;

		void countHeight(int& height) override;

		void rangeQuery(const Envelope& rect, std::vector<Feature>& features) override;

		bool NNQuery(double x, double y, std::vector<Feature>& features) override;

		RNode* pointInLeafNode(double x, double y) {
			if (root != nullptr) return root->pointInLeafNode(x, y);
			else return nullptr;
		}

		// Spatial Join: 输出满足空间距离条件的所有几何特征对
		// 正确性证明：
		// 1. Filter阶段：对于集合A中的每个几何对象a，扩展其包围盒dist距离
		// 2. 使用R树索引在集合B中查找与扩展包围盒相交的候选对象
		// 3. Refine阶段：对候选集做精确距离判断，保证不遗漏任何满足条件的对
		// 4. 由于包围盒扩展了dist距离，所有真正距离<=dist的对象必在候选集中
		void spatialJoin(const std::vector<Feature>& A,
			const std::vector<Feature>& B,
			double dist,
			std::vector<std::pair<Feature, Feature>>& result);

		// [[optional]] k最邻近几何特征查询 (k-NN)
		// 返回距离查询点最近的k个几何特征
		bool kNNQuery(double x, double y, int k, std::vector<Feature>& features);

		void draw() override { if (root != nullptr) root->draw(); }

		const Envelope& getEnvelope() const { return bbox; }

	public:
		static void test(int t);
		static void analyse();
	};

} // namespace hw6

#endif // !RTREE_H_INCLUDED
