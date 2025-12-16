#include "QuadTree.h"
#include <set>
#include <unordered_set>
#include <utility>
#include <limits>

namespace hw6 {

	/*
	 * QuadNode
	 */
	void QuadNode::split(size_t capacity) {
		for (int i = 0; i < 4; ++i) {
			delete children[i];
			children[i] = nullptr;
		}

		// Task construction
		// TODO
		double halfWidth = bbox.getWidth() / 2.0;
		double halfHeight = bbox.getHeight() / 2.0;
		double midX = bbox.getMinX() + halfWidth;
		double midY = bbox.getMinY() + halfHeight;

		// 1. 创建四个子节点的边界 (bbox)
		// NW (Top-Left): MinX, MidX, MidY, MaxY
		Envelope NW_bbox(bbox.getMinX(), midX, midY, bbox.getMaxY());
		// NE (Top-Right): MidX, MaxX, MidY, MaxY
		Envelope NE_bbox(midX, bbox.getMaxX(), midY, bbox.getMaxY());
		// SW (Bottom-Left): MinX, MidX, MinY, MidY
		Envelope SW_bbox(bbox.getMinX(), midX, bbox.getMinY(), midY);
		// SE (Bottom-Right): MidX, MaxX, MinY, MidY
		Envelope SE_bbox(midX, bbox.getMaxX(), bbox.getMinY(), midY);

		// 2. 初始化四个子节点 (它们都是叶子节点)
		children[0] = new QuadNode(NW_bbox);
		children[1] = new QuadNode(NE_bbox);
		children[2] = new QuadNode(SW_bbox);
		children[3] = new QuadNode(SE_bbox);

		// 3. 将当前节点中的所有Feature重新分配给子节点
		// 改进策略：
		// - 仅当某个子节点“完全包含”Feature的MBR时将其移动到该子节点
		// - 如果没有任何子节点完全包含该Feature，或被多个子节点包含，则保留在父节点（避免丢失或重复）
		std::vector<Feature> remaining;
		remaining.reserve(features.size());

		for (const Feature& feature : features) {
			Envelope featureBBox = feature.getEnvelope();

			int containCount = 0;
			int lastIndex = -1;
			for (int i = 0; i < 4; ++i) {
				if (children[i]->getEnvelope().contain(featureBBox)) {
					++containCount;
					lastIndex = i;
				}
			}

			if (containCount == 1) {
				// 唯一被某个子节点完全包含，移动到该子节点
				children[lastIndex]->add(feature);
			}
			else {
				// 无子节点或多个子节点包含 -> 保留在父节点
				remaining.push_back(feature);
			}
		}

		// 4. 将无法分配的要素保留在父节点
		features = std::move(remaining);

		// 5. 如果某个子节点的要素超过 capacity，可以递归分裂（使用传入的 capacity）
		for (int i = 0; i < 4; ++i) {
			if (children[i] && children[i]->features.size() > capacity) {
				children[i]->split(capacity);
			}
		}
		// isLeafNode() 将返回 false
	
	}

	void hw6::QuadNode::insert(const Feature& feature, size_t capacity) {
		const Envelope& fe = feature.getEnvelope();
		if (!bbox.contain(fe))
			return;

		if (isLeafNode()) {
			features.push_back(feature);
			// 使用传入的 capacity 参数
			if (capacity > 0 && features.size() > capacity) {
				split(capacity);
			}
			return;
		}

		int containCount = 0;
		int lastIndex = -1;
		for (int i = 0; i < 4; ++i) {
			if (children[i]->getEnvelope().contain(fe)) {
				++containCount;
				lastIndex = i;
			}
		}

		if (containCount == 1) {
			children[lastIndex]->insert(feature, capacity);
		} else {
			features.push_back(feature);
		}
	}

	void QuadNode::countNode(int& interiorNum, int& leafNum) {
		if (isLeafNode()) {
			++leafNum;
		}
		else {
			++interiorNum;
			for (int i = 0; i < 4; ++i)
				children[i]->countNode(interiorNum, leafNum);
		}
	}

	int QuadNode::countHeight(int height) {
		++height;
		if (!isLeafNode()) {
			int cur = height;
			for (int i = 0; i < 4; ++i) {
				height = std::max(height, children[i]->countHeight(cur));
			}
		}
		return height;
	}

	void QuadNode::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		if (!bbox.intersect(rect))
			return;

		// Task range query
		// TODO
		if (isLeafNode()) {
			// 1. 如果是叶子节点，遍历其所有Feature
			for (const Feature& feature : this->features) {
				// 2. 粗粒度过滤：检查Feature的MBR是否与查询矩形相交
				if (feature.getEnvelope().intersect(rect)) {
					// 3. 收集到候选集
					features.push_back(feature);
				}
			}
		}
		else {
			// 4. 如果是内部节点，递归调用子节点的rangeQuery
			for (int i = 0; i < 4; ++i) {
				children[i]->rangeQuery(rect, features);
			}
		}
	}

	QuadNode* QuadNode::pointInLeafNode(double x, double y) {
		// Task NN query
		// TODO
		// 1. 剪枝：如果当前节点不包含该点，返回空
		if (!bbox.contain(x, y)) {
			return nullptr;
		}

		// 2. 如果是叶子节点，返回自身
		if (isLeafNode()) {
			return this;
		}

		// 3. 如果是内部节点，递归查找子节点
		for (int i = 0; i < 4; ++i) {
			QuadNode* leaf = children[i]->pointInLeafNode(x, y);
			if (leaf) {
				return leaf;
			}
		}
		return nullptr;		//保留更安全
	}

	void QuadNode::draw() {
		if (isLeafNode()) {
			bbox.draw();
		}
		else {
			for (int i = 0; i < 4; ++i)
				children[i]->draw();
		}
	}

	/*
	 * QuadTree
	 */
	bool QuadTree::constructTree(const std::vector<Feature>& features) {
		if (features.empty())
			return false;

		// Task construction
		// TODO
		// 1. 计算所有 features 的联合包围盒
		Envelope unionBBox = features[0].getEnvelope();
		for (size_t i = 1; i < features.size(); ++i) {
			unionBBox = unionBBox.unionEnvelope(features[i].getEnvelope());
		}

		// 2. 确定根节点（QuadTree）的最终包围盒 (通常会扩大到方形)
		double minX = unionBBox.getMinX();
		double maxX = unionBBox.getMaxX();
		double minY = unionBBox.getMinY();
		double maxY = unionBBox.getMaxY();

		double width = maxX - minX;
		double height = maxY - minY;
		double maxDim = std::max(width, height);

		// 扩大到正方形以避免长条形MBR导致的低效划分
		// 居中扩展
		double expandX = (maxDim - width) / 2.0;
		double expandY = (maxDim - height) / 2.0;

		bbox = Envelope(minX - expandX, maxX + expandX, minY - expandY, maxY + expandY);

		// 3. 初始化根节点
		if (root) {
			delete root;
			root = nullptr;
		}
		// 使用计算出的、可能扩大后的包围盒作为根节点范围
		root = new QuadNode(bbox, capacity);

		// 4. 将所有 features 插入到根节点
		for (const Feature& feature : features) {
			root->insert(feature,capacity);
		}

		// 兼容性：如果 insert 未触发任何分裂（根仍为叶子）且输入数量超过 capacity，则强制 split 一次
		if (features.size() > capacity && root->isLeafNode()) {
			root->split(capacity);
		}

		//bbox = Envelope(-74.1, -73.8, 40.6, 40.8); // 注意此行代码需要更新为features的包围盒，或根节点的包围盒

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

		// Task range query
		// TODO

		// filter step (选择查询区域与几何对象包围盒相交的几何对象)
		if (root) {
			root->rangeQuery(rect, features);
		}
		// 注意四叉树区域查询仅返回候选集，精炼步在hw6的rangeQuery中完成
	}

	bool QuadTree::NNQuery(double x, double y, std::vector<Feature>& features) {
		if (!root || !(root->getEnvelope().contain(x, y)))
			return false;

		// Task NN query
		// TODO

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

		// 注意：这里保持与 hw6.cpp 的约定——索引只返回候选集，精炼与去重留给 hw6.cpp
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
