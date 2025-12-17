#ifndef QUADTREE_H_INCLUDED
#define QUADTREE_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"
#include <string>

namespace hw6 {

	class QuadNode {
	private:
		Envelope bbox;
		QuadNode* children[4];
		std::vector<Feature> features;
		size_t nodeCapacity;

	public:
		QuadNode() = delete;
		QuadNode(const Envelope& box) : bbox(box), nodeCapacity(0) {
			children[0] = children[1] = children[2] = children[3] = nullptr;
		}

		QuadNode(const Envelope& box, size_t capacity) : bbox(box), nodeCapacity(capacity) {
			children[0] = children[1] = children[2] = children[3] = nullptr;
		}

		~QuadNode() {
			for (int i = 0; i < 4; ++i) {
				delete children[i];
				children[i] = nullptr;
			}
		}

		bool isLeafNode() { return children[0] == nullptr; }

		const Envelope& getEnvelope() { return bbox; }

		QuadNode* getChildNode(size_t i) { return i < 4 ? children[i] : nullptr; }

		size_t getFeatureNum() const { return features.size(); }

		const Feature& getFeature(size_t i) const { return features[i]; }

		void add(const Feature& f) { features.push_back(f); }

		void add(const std::vector<Feature>& fs) {
			features.insert(features.begin(), fs.begin(), fs.end());
		}

		// 插入接口（使用传入 capacity）
		void insert(const Feature& feature, size_t capacity);

		void countNode(int& interiorNum, int& leafNum);

		int countHeight(int height);

		void draw();

		// split the node into four child nodes, assign each feature to its
		// overlaped child node(s), clear feature vector, and split child node(s) if
		// its number of features is larger than capacity
		void split(size_t capacity);

		void rangeQuery(const Envelope& rect, std::vector<Feature>& features);

		QuadNode* pointInLeafNode(double x, double y);
	};

	class QuadTree : public Tree {
	private:
		QuadNode* root;

	public:
		QuadTree() : Tree(5), root(nullptr) {}
		QuadTree(size_t cap) : Tree(cap), root(nullptr) {}
		~QuadTree() {
			if (root != nullptr)
				delete root;
			root = nullptr;
		}

		virtual bool constructTree(const std::vector<Feature>& features) override;

		virtual void countNode(int& interiorNum, int& leafNum) override;

		virtual void countHeight(int& height) override;

		virtual void rangeQuery(const Envelope& rect, std::vector<Feature>& features) override;

		virtual bool NNQuery(double x, double y, std::vector<Feature>& features) override;

		QuadNode* pointInLeafNode(double x, double y) {
			return root->pointInLeafNode(x, y);
		}

		virtual void draw() override;
		void spatialJoin(const std::vector<Feature>& A,
			const std::vector<Feature>& B,
			double dist,
			std::vector<std::pair<Feature, Feature>>& result);

	public:
		static void test(int t);

		static void analyse();
	};



} // namespace hw6

#endif