#include "QuadTree.h"
#include <set>

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
		children[NW] = new QuadNode(NW_bbox);
		children[NE] = new QuadNode(NE_bbox);
		children[SW] = new QuadNode(SW_bbox);
		children[SE] = new QuadNode(SE_bbox);

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
				children[lastIndex]->addFeature(feature);
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

	// 增加：实现 insert 方法
	void QuadNode::insert(const Feature& feature) {
		// 如果 feature 的 MBR 不在当前节点 bbox 内，直接返回（上层应保证 root 覆盖所有 features）
		const Envelope& fe = feature.getEnvelope();
		if (!bbox.contain(fe))
			return;

		// 若为叶子节点，尝试加入并根据 nodeCapacity 分裂
		if (isLeafNode()) {
			features.push_back(feature);
			// 使用节点存储的 nodeCapacity，如果为 0 则不触发自动分裂（保持向后兼容）
			if (nodeCapacity > 0 && features.size() > nodeCapacity) {
				split(nodeCapacity);
			}
			return;
		}

		// 非叶节点，尝试下推到唯一完全包含的子节点，否则保留在当前节点
		int containCount = 0;
		int lastIndex = -1;
		for (int i = 0; i < 4; ++i) {
			if (children[i]->getEnvelope().contain(fe)) {
				++containCount;
				lastIndex = i;
			}
		}

		if (containCount == 1) {
			children[lastIndex]->insert(feature);
		}
		else {
			// 无子节点或多个子节点包含 -> 存在于父节点
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
			root->insert(feature);
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

		// 1. 遍历四叉树所有节点，计算 min over features of maxDistance2Envelope(x,y)
		double minMaxDist = std::numeric_limits<double>::infinity();

		std::function<void(QuadNode*)> visit = [&](QuadNode* node) {
			if (!node) return;

			// 处理当前节点存储的 features（父节点也可能保留若干 feature）
			size_t fn = node->getFeatureNum();
			for (size_t i = 0; i < fn; ++i) {
				const Feature& f = node->getFeature(i);
				double d = f.maxDistance2Envelope(x, y);
				if (d < minMaxDist) minMaxDist = d;
			}

			// 递归子节点
			if (!node->isLeafNode()) {
				for (size_t i = 0; i < 4; ++i) {
					QuadNode* c = node->getChildNode(i);
					if (c) visit(c);
				}
			}
			};

		visit(root);

		if (!std::isfinite(minMaxDist))
			return false;

		// 2. 用得到的最短最大距离构造查询矩形并用 rangeQuery 得到候选集
		Envelope qbox(x - minMaxDist, x + minMaxDist, y - minMaxDist, y + minMaxDist);
		rangeQuery(qbox, features);

		// 返回是否有候选；精确最近邻由上层 hw6::NNQuery 负责
		return !features.empty();
	}

	void QuadTree::draw() {
		if (root)
			root->draw();
	}

} // namespace hw6
