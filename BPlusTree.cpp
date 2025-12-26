#include "BPlusTree.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>
#include <unordered_set>

namespace hw6 {

	/*
	 * SpaceFilling 实现
	 */
	uint64_t SpaceFilling::xyToZOrder(uint32_t x, uint32_t y) {
		// Morton编码：交错x和y的比特位
		uint64_t z = 0;
		for (int i = 0; i < 32; ++i) {
			z |= ((uint64_t)(x & (1u << i)) << i) | ((uint64_t)(y & (1u << i)) << (i + 1));
		}
		return z;
	}

	uint32_t SpaceFilling::normalizeCoord(double coord, double minCoord, double maxCoord, uint32_t gridSize) {
		if (maxCoord <= minCoord) return 0;
		double normalized = (coord - minCoord) / (maxCoord - minCoord);
		normalized = std::max(0.0, std::min(1.0, normalized));
		return static_cast<uint32_t>(normalized * (gridSize - 1));
	}

	uint64_t SpaceFilling::computeZValue(const Feature& feature, const Envelope& bbox, uint32_t gridSize) {
		const Envelope& fEnv = feature.getEnvelope();
		// 使用包围盒中心点计算Z值
		double centerX = (fEnv.getMinX() + fEnv.getMaxX()) / 2.0;
		double centerY = (fEnv.getMinY() + fEnv.getMaxY()) / 2.0;

		uint32_t gridX = normalizeCoord(centerX, bbox.getMinX(), bbox.getMaxX(), gridSize);
		uint32_t gridY = normalizeCoord(centerY, bbox.getMinY(), bbox.getMaxY(), gridSize);

		return xyToZOrder(gridX, gridY);
	}

	/*
	 * BPlusNode 实现
	 */
	BPlusNode::BPlusNode(bool leaf) : isLeaf(leaf), parent(nullptr), next(nullptr) {}

	BPlusNode::~BPlusNode() {
		for (auto* child : children) {
			delete child;
		}
		children.clear();
	}

	size_t BPlusNode::findChildIndex(uint64_t key) const {
		// 二分查找第一个 >= key 的位置
		size_t left = 0, right = keys.size();
		while (left < right) {
			size_t mid = left + (right - left) / 2;
			if (keys[mid] < key) {
				left = mid + 1;
			}
			else {
				right = mid;
			}
		}
		return left;
	}

	void BPlusNode::insertLeaf(uint64_t key, const Feature& feature, size_t capacity, BPlusNode*& newNode) {
		newNode = nullptr;

		// 找到插入位置
		size_t pos = findChildIndex(key);

		// 插入键值对
		keys.insert(keys.begin() + pos, key);
		features.insert(features.begin() + pos, feature);

		// 检查是否需要分裂
		if (keys.size() > capacity) {
			split(capacity, newNode, key);
		}
	}

	void BPlusNode::insertInternal(uint64_t key, BPlusNode* child) {
		size_t pos = findChildIndex(key);
		keys.insert(keys.begin() + pos, key);
		children.insert(children.begin() + pos + 1, child);
		child->setParent(this);
	}

	void BPlusNode::split(size_t capacity, BPlusNode*& newNode, uint64_t& midKey) {
		size_t mid = (keys.size() + 1) / 2;
		midKey = keys[mid];

		newNode = new BPlusNode(isLeaf);

		if (isLeaf) {
			// 叶节点分裂：右半部分移动到新节点
			newNode->keys.assign(keys.begin() + mid, keys.end());
			newNode->features.assign(features.begin() + mid, features.end());
			
			keys.erase(keys.begin() + mid, keys.end());
			features.erase(features.begin() + mid, features.end());

			// 更新叶节点链表
			newNode->next = this->next;
			this->next = newNode;
		}
		else {
			// 内部节点分裂
			newNode->keys.assign(keys.begin() + mid + 1, keys.end());
			newNode->children.assign(children.begin() + mid + 1, children.end());

			keys.erase(keys.begin() + mid, keys.end());
			children.erase(children.begin() + mid + 1, children.end());

			// 更新父节点指针
			for (auto* child : newNode->children) {
				child->setParent(newNode);
			}
		}
	}

	/*
	 * BPlusTree 实现
	 */
	BPlusTree::BPlusTree() : Tree(20), root(nullptr), gridSize(1024) {}

	BPlusTree::BPlusTree(size_t cap, uint32_t grid) : Tree(cap), root(nullptr), gridSize(grid) {}

	BPlusTree::~BPlusTree() {
		delete root;
		root = nullptr;
	}

	bool BPlusTree::constructTree(const std::vector<Feature>& features) {
		if (features.empty()) return false;

		// 1. 计算全局包围盒
		bbox = features[0].getEnvelope();
		for (size_t i = 1; i < features.size(); ++i) {
			bbox = bbox.unionEnvelope(features[i].getEnvelope());
		}

		// 2. 计算每个特征的Z值并排序
		std::vector<std::pair<uint64_t, Feature>> zFeatures;
		zFeatures.reserve(features.size());

		for (const Feature& f : features) {
			uint64_t z = SpaceFilling::computeZValue(f, bbox, gridSize);
			zFeatures.emplace_back(z, f);
		}

		// 按Z值排序
		std::sort(zFeatures.begin(), zFeatures.end(),
			[](const std::pair<uint64_t, Feature>& a, const std::pair<uint64_t, Feature>& b) {
				return a.first < b.first;
			});

		// 3. 构建B+树
		delete root;
		root = new BPlusNode(true);

		for (const auto& zf : zFeatures) {
			BPlusNode* newNode = nullptr;
			uint64_t midKey = 0;
			insertHelper(zf.first, zf.second, root, newNode, midKey);

			// 如果根节点分裂，创建新根
			if (newNode != nullptr) {
				BPlusNode* newRoot = new BPlusNode(false);
				// 使用新增的公共方法
				newRoot->addKey(midKey);
				newRoot->addChild(root);
				newRoot->addChild(newNode);
				root = newRoot;
			}
		}

		return true;
	}

	void BPlusTree::insertHelper(uint64_t zValue, const Feature& feature, BPlusNode* node, BPlusNode*& newNode, uint64_t& midKey) {
		if (node->getIsLeaf()) {
			// 叶节点：直接插入
			node->insertLeaf(zValue, feature, capacity, newNode);
			if (newNode) {
				midKey = newNode->getKey(0);
			}
		}
		else {
			// 内部节点：递归插入
			size_t childIndex = node->findChildIndex(zValue);
			BPlusNode* child = node->getChild(childIndex);

			BPlusNode* newChild = nullptr;
			uint64_t childMidKey = 0;
			insertHelper(zValue, feature, child, newChild, childMidKey);

			// 如果子节点分裂，插入新键到当前节点
			if (newChild != nullptr) {
				node->insertInternal(childMidKey, newChild);

				// 检查当前节点是否需要分裂
				if (node->getKeyNum() > capacity) {
					node->split(capacity, newNode, midKey);
				}
			}
		}
	}

	void BPlusTree::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		features.clear();
		if (!root) return;

		// 计算查询区域的Z值范围
		// 由于Z-order保持空间局部性，我们遍历叶节点链表
		BPlusNode* leaf = root;
		while (!leaf->getIsLeaf()) {
			leaf = leaf->getChild(0);
		}

		// 遍历所有叶节点
		while (leaf) {
			for (size_t i = 0; i < leaf->getFeatureNum(); ++i) {
				const Feature& f = leaf->getFeature(i);
				if (f.getEnvelope().intersect(rect)) {
					features.push_back(f);
				}
			}
			leaf = leaf->getNext();
		}
	}

	bool BPlusTree::NNQuery(double x, double y, std::vector<Feature>& features) {
		features.clear();
		if (!root) return false;

		// 计算查询点的Z值
		Envelope pointEnv(x, x, y, y);
		Point queryPoint(x, y);
		uint64_t queryZ = SpaceFilling::computeZValue(Feature("", &queryPoint), bbox, gridSize);

		// 找到包含查询点Z值的叶节点
		BPlusNode* targetLeaf = root;
		while (!targetLeaf->getIsLeaf()) {
			size_t childIndex = targetLeaf->findChildIndex(queryZ);
			targetLeaf = targetLeaf->getChild(childIndex);
		}

		// 计算最小的maxDistance2Envelope
		double minMaxDist = std::numeric_limits<double>::infinity();

		// 从目标叶节点开始，向两侧扩展搜索
		BPlusNode* leaf = root;
		while (!leaf->getIsLeaf()) {
			leaf = leaf->getChild(0);
		}

		while (leaf) {
			for (size_t i = 0; i < leaf->getFeatureNum(); ++i) {
				const Feature& f = leaf->getFeature(i);
				double d = f.maxDistance2Envelope(x, y);
				if (d < minMaxDist) {
					minMaxDist = d;
				}
			}
			leaf = leaf->getNext();
		}

		if (!std::isfinite(minMaxDist)) {
			return false;
		}

		// 使用区域查询获取候选集
		Envelope queryRect(x - minMaxDist, x + minMaxDist, y - minMaxDist, y + minMaxDist);
		rangeQuery(queryRect, features);

		return !features.empty();
	}

	void BPlusTree::spatialJoin(const std::vector<Feature>& A,
		const std::vector<Feature>& B,
		double dist,
		std::vector<std::pair<Feature, Feature>>& result) {

		result.clear();
		if (A.empty() || B.empty()) return;

		// 构建集合B的B+树索引
		BPlusTree treeB(this->capacity, this->gridSize);
		if (!treeB.constructTree(B)) {
			return;
		}

		std::vector<Feature> candidates;
		candidates.reserve(64);

		for (const Feature& a : A) {
			const Envelope& ea = a.getEnvelope();
			Envelope searchEnv(ea.getMinX() - dist, ea.getMaxX() + dist,
				ea.getMinY() - dist, ea.getMaxY() + dist);

			// Filter: 获取候选集
			treeB.rangeQuery(searchEnv, candidates);

			// Refine: 精确距离计算和去重
			std::unordered_set<const Geometry*> seenB;
			seenB.reserve(candidates.size() * 2);

			for (const Feature& b : candidates) {
				const Geometry* ga = a.getGeom();
				const Geometry* gb = b.getGeom();
				if (!ga || !gb) continue;

				if (!seenB.insert(gb).second) continue;

				double d = ga->distance(gb);
				if (d <= dist) {
					result.emplace_back(a, b);
				}
			}
		}
	}

	void BPlusTree::countNode(int& interiorNum, int& leafNum) {
		interiorNum = 0;
		leafNum = 0;
		if (root) {
			countNodeHelper(root, interiorNum, leafNum);
		}
	}

	void BPlusTree::countNodeHelper(BPlusNode* node, int& interiorNum, int& leafNum) {
		if (!node) return;

		if (node->getIsLeaf()) {
			++leafNum;
		}
		else {
			++interiorNum;
			for (size_t i = 0; i < node->getChildNum(); ++i) {
				countNodeHelper(node->getChild(i), interiorNum, leafNum);
			}
		}
	}

	void BPlusTree::countHeight(int& height) {
		height = 0;
		if (root) {
			height = countHeightHelper(root);
		}
	}

	int BPlusTree::countHeightHelper(BPlusNode* node) {
		if (!node || node->getIsLeaf()) {
			return 1;
		}

		int maxHeight = 0;
		for (size_t i = 0; i < node->getChildNum(); ++i) {
			int h = countHeightHelper(node->getChild(i));
			maxHeight = std::max(maxHeight, h);
		}
		return maxHeight + 1;
	}

	void BPlusTree::draw() {
		// 可选：实现绘制功能
		// 这里留空，因为B+树的可视化较为复杂
	}


} // namespace hw6