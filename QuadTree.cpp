#include "QuadTree.h"
#include <set>
#include <unordered_set>
#include <utility>
#include <limits>

namespace hw6 {

	/*
	 * QuadNode
	 */
	void QuadNode::insert(const Feature& feature, size_t capacity) {
		if (isLeafNode()) {
			// 叶节点:直接添加feature
			add(feature);
			// 检查是否超容量，需要分裂
			if (features.size() > capacity) {
				split(capacity);
			}						
		}
		else {
			// 内部节点:找到与feature包围盒相交的子节点并递归插入
			const Envelope& fEnv = feature.getEnvelope();
			for (int i = 0; i < 4; ++i) {
				if (children[i] && children[i]->getEnvelope().intersect(fEnv)) {
					children[i]->insert(feature, capacity);
				}
			}
		}
	}

	void QuadNode::split(size_t capacity) {
		// 清空现有子节点
		for (int i = 0; i < 4; ++i) {
			delete children[i];
			children[i] = nullptr;
		}

		// 计算四个象限的包围盒
		double minX = bbox.getMinX();
		double maxX = bbox.getMaxX();
		double minY = bbox.getMinY();
		double maxY = bbox.getMaxY();
		double midX = (minX + maxX) / 2.0;
		double midY = (minY + maxY) / 2.0;

		// 创建四个子节点 (西南、东南、西北、东北)
		children[0] = new QuadNode(Envelope(minX, midX, minY, midY), capacity); // 西南
		children[1] = new QuadNode(Envelope(midX, maxX, minY, midY), capacity); // 东南
		children[2] = new QuadNode(Envelope(minX, midX, midY, maxY), capacity); // 西北
		children[3] = new QuadNode(Envelope(midX, maxX, midY, maxY), capacity); // 东北

		// 将当前节点的features分配到相应的子节点
		for (const Feature& f : features) {
			const Envelope& fEnv = f.getEnvelope();
			for (int i = 0; i < 4; ++i) {
				if (children[i]->getEnvelope().intersect(fEnv)) {
					children[i]->add(f);
				}
			}
		}

		// 清空当前节点的features
		features.clear();

		// 递归分裂超容量的子节点
		for (int i = 0; i < 4; ++i) {
			if (children[i]->getFeatureNum() > capacity) {
				children[i]->split(capacity);
			}
		}
	}

	void QuadNode::countNode(int& interiorNum, int& leafNum) {
		if (isLeafNode()) {
			++leafNum;
		}
		else {
			++interiorNum;
			for (int i = 0; i < 4; ++i) {
				if (children[i]) {
					children[i]->countNode(interiorNum, leafNum);
				}
			}
		}
	}

	int QuadNode::countHeight(int height) {
		++height;
		if (!isLeafNode()) {
			int maxHeight = height;
			for (int i = 0; i < 4; ++i) {
				if (children[i]) {
					maxHeight = std::max(maxHeight, children[i]->countHeight(height));
				}
			}
			return maxHeight;
		}
		return height;
	}

	void QuadNode::draw() {
		if (isLeafNode()) {
			bbox.draw();
		}
		else {
			for (int i = 0; i < 4; ++i) {
				if (children[i]) {
					children[i]->draw();
				}
			}
		}
	}

	void QuadNode::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		if (!bbox.intersect(rect))
			return;

		if (isLeafNode()) {
			// 叶节点:检查每个feature的包围盒是否与查询区域相交
			for (const Feature& f : this->features) {
				if (f.getEnvelope().intersect(rect)) {
					features.push_back(f);
				}
			}
		}
		else {
			// 内部节点:递归查询相交的子节点
			for (int i = 0; i < 4; ++i) {
				if (children[i]) {
					children[i]->rangeQuery(rect, features);
				}
			}
		}
	}

	QuadNode* QuadNode::pointInLeafNode(double x, double y) {
		if (!bbox.contain(x, y)) {
			return nullptr;
		}

		if (isLeafNode()) {
			return this;
		}
		else {
			for (int i = 0; i < 4; ++i) {
				if (children[i]) {
					QuadNode* result = children[i]->pointInLeafNode(x, y);
					if (result != nullptr) {
						return result;
					}
				}
			}
		}
		return nullptr;
	}

	/*
	 * QuadTree
	 */
	bool QuadTree::constructTree(const std::vector<Feature>& features) {
		if (features.empty())
			return false;

		// 1. 计算所有 features 的联合包围盒
		Envelope unionBBox = features[0].getEnvelope();
		for (size_t i = 1; i < features.size(); ++i) {
			unionBBox = unionBBox.unionEnvelope(features[i].getEnvelope());
		}

		// 2. 设置树的包围盒
		bbox = unionBBox;

		// 3. 初始化根节点
		if (root) {
			delete root;
			root = nullptr;
		}
		root = new QuadNode(bbox, capacity);

		// 4. 逐个插入 features
		for (const Feature& feature : features) {
			root->insert(feature, capacity);
		}

		return true;
	}

	void QuadTree::countNode(int& interiorNum, int& leafNum) {
		interiorNum = 0;
		leafNum = 0;
		if (root)
			root->countNode(interiorNum, leafNum);
	}

	void QuadTree::countHeight(int& height) {
		height = 0;
		if (root)
			height = root->countHeight(0);
	}

	void QuadTree::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		features.clear();

		if (root) {
			root->rangeQuery(rect, features);
		}
	}

	bool QuadTree::NNQuery(double x, double y, std::vector<Feature>& features) {
		if (!root || !(root->getEnvelope().contain(x, y)))
			return false;

		features.clear();

		// 计算所有节点的最小 maxDistance2Envelope（遍历整棵树）
		double minMaxDist = std::numeric_limits<double>::infinity();

		// 使用显式栈遍历所有节点
		std::vector<QuadNode*> stack;
		stack.reserve(64);
		stack.push_back(root);

		while (!stack.empty()) {
			QuadNode* node = stack.back();
			stack.pop_back();
			if (!node) continue;

			for (size_t i = 0; i < node->getFeatureNum(); ++i) {
				const Feature& f = node->getFeature(i);
				double d = f.maxDistance2Envelope(x, y);
				if (d < minMaxDist) minMaxDist = d;
			}

			if (!node->isLeafNode()) {
				for (size_t i = 0; i < 4; ++i) {
					QuadNode* c = node->getChildNode(i);
					if (c) stack.push_back(c);
				}
			}
		}

		if (!std::isfinite(minMaxDist)) {
			features.clear();
			return false;
		}

		Envelope qbox(x - minMaxDist, x + minMaxDist, y - minMaxDist, y + minMaxDist);
		// 调用已有的 rangeQuery（仅做 MBR 粗筛），返回候选集到 features
		this->rangeQuery(qbox, features);

		// 索引只返回候选集，精炼与去重留给 hw6.cpp
		return !features.empty();
	}

	void QuadTree::draw() {
		if (root)
			root->draw();
	}

	void QuadTree::spatialJoin(const std::vector<Feature>& A,
		const std::vector<Feature>& B,
		double dist,
		std::vector<std::pair<Feature, Feature>>& result) {

		result.clear();

		if (A.empty() || B.empty())
			return;

		// 构建用于 B 的四叉树索引（局部临时）
		QuadTree treeB(this->capacity);
		treeB.setCapacity(this->capacity);
		if (!treeB.constructTree(B)) {
			// 若构建失败（例如 B 为空），直接返回
			return;
		}

		std::vector<Feature> candidates;
		candidates.reserve(64);

		for (const Feature& a : A) {
			candidates.clear(); // 每个 a 单独收集候选，避免累加

			// 扩展 a 的包围盒
			const Envelope& eb = a.getEnvelope();
			Envelope searchEnv(eb.getMinX() - dist, eb.getMaxX() + dist,
				eb.getMinY() - dist, eb.getMaxY() + dist);

			// filter: 用 treeB 获取与扩展矩形相交的候选集
			treeB.rangeQuery(searchEnv, candidates);

			// refine: 对候选集做精确距离判断，并去重（基于 Geometry 指针）
			std::unordered_set<const Geometry*> seenB;
			seenB.reserve(candidates.size() * 2);

			for (const Feature& b : candidates) {
				const Geometry* gb = b.getGeom();
				const Geometry* ga = a.getGeom();
				if (!ga || !gb) continue;

				// 去重同一几何
				if (!seenB.insert(gb).second) continue;

				double d = ga->distance(gb);
				if (d <= dist) {
					result.emplace_back(a, b);
				}
			}
		}
	}

} // namespace hw6
