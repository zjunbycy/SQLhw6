#include "RTree.h"
#include <random>

namespace hw6 {

	// RNode 实现
	void RNode::add(RNode* child) {
		children.push_back(child);
		child->parent = this;
		++childrenNum;
	}

	void RNode::remove(const Feature& f) {
		auto where = [&]() {
			for (auto itr = features.begin(); itr != features.end(); ++itr)
				if (itr->getName() == f.getName())
					return itr;
			return features.end();
			}();
		if (where != features.end()) {
			features.erase(where);
		}
		if (features.empty())
			features.shrink_to_fit(); // free memory unused but allocated
	}

	void RNode::remove(RNode* child) {
		for (int i = 0; i < childrenNum; ++i)
			if (children[i] == child) {
				--childrenNum;
				std::swap(children[i], children[childrenNum]);
				children.pop_back();
				break;
			}
	}

	Feature RNode::popBackFeature() {
		auto ret = features.back();
		features.pop_back();
		return ret;
	}

	RNode* RNode::popBackChildNode() {
		if (childrenNum == 0) return nullptr;
		--childrenNum;
		auto ret = children[childrenNum];
		children.pop_back();
		return ret;
	}

	void RNode::countNode(int& interiorNum, int& leafNum) {
		if (isLeafNode()) {
			++leafNum;
		}
		else {
			++interiorNum;
			for (int i = 0; i < childrenNum; ++i)
				children[i]->countNode(interiorNum, leafNum);
		}
	}

	int RNode::countHeight(int height) {
		++height;
		if (!isLeafNode()) {
			int cur = height;
			for (int i = 0; i < childrenNum; ++i)
				height = std::max(height, children[i]->countHeight(cur));
		}
		return height;
	}

	void RNode::draw() {
		if (isLeafNode()) {
			bbox.draw();
		}
		else
			for (int i = 0; i < childrenNum; ++i)
				children[i]->draw();
	}

	void RNode::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		// Task rangeQuery
		/* TODO */
		if (!bbox.intersect(rect))
			return;

		if (isLeafNode()) {
			// 叶节点：检查每个Feature的包围盒是否与查询区域相交
			for (size_t i = 0; i < this->features.size(); ++i) {
				const Feature& f = this->features[i];
				if (f.getEnvelope().intersect(rect)) {
					features.push_back(f);
				}
			}
		}
		else {
			// 内部节点：递归查询每个子节点
			for (int i = 0; i < childrenNum; ++i)
				children[i]->rangeQuery(rect, features);
		}
		// filter step (选择查询区域与几何对象包围盒相交的几何对象)
		// 注意R树区域查询仅返回候选集，精炼步在hw6的rangeQuery中完成
	}

	RNode* RNode::pointInLeafNode(double x, double y) {
		// Task pointInLeafNode
		/* TODO */
		if(!bbox.contain(x,y)){return nullptr;}
		if(isLeafNode()){
			return this;
		}else{
			for(int i=0;i<childrenNum;++i){
				RNode* res=children[i]->pointInLeafNode(x,y);
				if(res!=nullptr){
					return res;
				}
			}
		}
		return nullptr;
	}


	// RTree 实现
	RTree::RTree(int maxChildren) : Tree(maxChildren), maxChildren(maxChildren) {
		if (maxChildren < 4) throw std::invalid_argument("maxChildren must be >= 4");
	}

	void RTree::countNode(int& interiorNum, int& leafNum) {
		interiorNum = leafNum = 0;
		if (root != nullptr)
			root->countNode(interiorNum, leafNum);
	}

	void RTree::countHeight(int& height) {
		height = 0;
		if (root != nullptr)
			height = root->countHeight(height);
	}

bool RTree::constructTree(const std::vector<Feature>& features) {
		// Task RTree construction
		/* TODO */
		// 几何特征按 x 轴顺序逐个插入；节点选择面积增量最小；
		// 溢出时采用二次分裂，种子为最左/最右；为增强平衡性，在相同 minX 的局部段内随机插入。
		//输入检查
		if (features.empty())
			return false;
		//面积增量
		auto areaEnlargement = [](const Envelope& nodeBox, const Envelope& box) {
			Envelope u = nodeBox.unionEnvelope(box);
			return u.getArea() - nodeBox.getArea();
		};
		//从当前内部节点选择子树（面积增量最小，若相同选择面积更小）
		auto chooseSubtreeLocal = [&](RNode* node, const Envelope& box) -> RNode* {
			RNode* best = nullptr;
			double minEnl = std::numeric_limits<double>::max();
			double minArea = std::numeric_limits<double>::max();
			for (int i = 0; i < node->getChildNum(); ++i) {
				RNode* c = node->getChildNode(static_cast<size_t>(i));
				const Envelope& cb = c->getEnvelope();
				double enl = areaEnlargement(cb, box);
				double area = cb.getArea();
				if (enl < minEnl || (enl == minEnl && area < minArea)) {minEnl = enl;minArea = area;best = c;}
			}return best;};
		//向上更新bbox
		auto updateBboxUp = [&](RNode* node) {
			while (node) {
				if (node->isLeafNode()) {
					if (node->getFeatureNum() > 0) {
						Envelope nb = node->getFeature(0).getEnvelope();
						for (size_t i = 1; i < node->getFeatureNum(); ++i)
							nb = nb.unionEnvelope(node->getFeature(i).getEnvelope());
						node->setEnvelope(nb);
					}
				} else {
					if (node->getChildNum() > 0) {
						Envelope nb = node->getChildNode(0)->getEnvelope();
						for (int i = 1; i < node->getChildNum(); ++i)
							nb = nb.unionEnvelope(node->getChildNode(static_cast<size_t>(i))->getEnvelope());
						node->setEnvelope(nb);
					}
				}
				node = node->getParent();
			}
		};

		// 叶节点二次分裂，种子为最左(minX最小)和最右(maxX最大)
		auto splitLeafQuadratic = [&](RNode* leaf, const Feature& newF, RNode*& newLeaf) {
			std::vector<Feature> all = leaf->getFeatures();
			all.push_back(newF);

			// 选择最左和最右作为种子
			int leftIdx = 0, rightIdx = 0;
			double minX = std::numeric_limits<double>::max();
			double maxX = -std::numeric_limits<double>::max();
			for (int i = 0; i < static_cast<int>(all.size()); ++i) {
				const Envelope& e = all[i].getEnvelope();
				if (e.getMinX() < minX) { minX = e.getMinX(); leftIdx = i; }
				if (e.getMaxX() > maxX) { maxX = e.getMaxX(); rightIdx = i; }
			}
			if (leftIdx == rightIdx) {
				// 退化情况：选择不同的另一个点作为右种子
				rightIdx = (rightIdx + 1) % static_cast<int>(all.size());
			}

			std::vector<Feature> g1, g2;
			g1.push_back(all[leftIdx]);g2.push_back(all[rightIdx]);
			Envelope b1 = all[leftIdx].getEnvelope();Envelope b2 = all[rightIdx].getEnvelope();

			std::vector<char> used(all.size(), 0);
			used[leftIdx] = used[rightIdx] = 1;

			// 使用 ceil(maxChildren/2.0) 作为最小条目数
			int minEntries = std::max(1, (maxChildren + 1) / 2);
			int assigned = 2;
			while (assigned < static_cast<int>(all.size())) {
				// 若某组必须接收剩余所有元素以满足最小数目
				int remaining = static_cast<int>(all.size()) - assigned;
				if (static_cast<int>(g1.size()) + remaining == minEntries) {
					for (int i = 0; i < static_cast<int>(all.size()); ++i)
						if (!used[i]) { g1.push_back(all[i]); b1 = b1.unionEnvelope(all[i].getEnvelope()); used[i] = 1; ++assigned; }
					break;
				}
				if (static_cast<int>(g2.size()) + remaining == minEntries) {
					for (int i = 0; i < static_cast<int>(all.size()); ++i)
						if (!used[i]) { g2.push_back(all[i]); b2 = b2.unionEnvelope(all[i].getEnvelope()); used[i] = 1; ++assigned; }
					break;
				}
				// 选择使两组面积增量差异最大的元素
				int choose = -1;
				double maxDiff = -std::numeric_limits<double>::max();
				for (int i = 0; i < static_cast<int>(all.size()); ++i) {
					if (used[i]) continue;
					const Envelope& eb = all[i].getEnvelope();
					double d1 = areaEnlargement(b1, eb);double d2 = areaEnlargement(b2, eb);double diff = std::abs(d1 - d2);
					if (diff > maxDiff) { maxDiff = diff; choose = i; }
				}
				if (choose == -1) break;
				const Envelope& eb = all[choose].getEnvelope();
				double d1 = areaEnlargement(b1, eb);double d2 = areaEnlargement(b2, eb);
				if (d1 < d2 || (d1 == d2 && g1.size() <= g2.size())) {
					g1.push_back(all[choose]); b1 = b1.unionEnvelope(eb);
				} else {
					g2.push_back(all[choose]); b2 = b2.unionEnvelope(eb);
				}
				used[choose] = 1; ++assigned;
			}

			// 重建当前叶节点和新叶节点
			while (leaf->getFeatureNum() > 0) leaf->popBackFeature();
			for (const auto& f : g1) leaf->add(f);
			leaf->setEnvelope(b1);

			newLeaf = new RNode(b2);
			for (const auto& f : g2) newLeaf->add(f);
			newLeaf->setEnvelope(b2);
		};

		// 内部节点二次分裂（对children分组，种子仍按最左/最右）
		auto splitInnerQuadratic = [&](RNode* inner, RNode* extraChild, RNode*& newInner) {
			std::vector<RNode*> all;
			for (int i = 0; i < inner->getChildNum(); ++i)
				all.push_back(inner->getChildNode(static_cast<size_t>(i)));
			all.push_back(extraChild);

			int leftIdx = 0, rightIdx = 0;
			double minX = std::numeric_limits<double>::max();
			double maxX = -std::numeric_limits<double>::max();
			for (int i = 0; i < static_cast<int>(all.size()); ++i) {
				const Envelope& e = all[i]->getEnvelope();
				if (e.getMinX() < minX) { minX = e.getMinX(); leftIdx = i; }
				if (e.getMaxX() > maxX) { maxX = e.getMaxX(); rightIdx = i; }
			}
			if (leftIdx == rightIdx) rightIdx = (rightIdx + 1) % static_cast<int>(all.size());

			std::vector<RNode*> g1, g2;
			g1.push_back(all[leftIdx]);g2.push_back(all[rightIdx]);
			Envelope b1 = all[leftIdx]->getEnvelope();Envelope b2 = all[rightIdx]->getEnvelope();

			std::vector<char> used(all.size(), 0);
			used[leftIdx] = used[rightIdx] = 1;

			int minEntries = std::max(1, (maxChildren + 1) / 2);
			int assigned = 2;
			while (assigned < static_cast<int>(all.size())) {
				int remaining = static_cast<int>(all.size()) - assigned;
				if (static_cast<int>(g1.size()) + remaining == minEntries) {
					for (int i = 0; i < static_cast<int>(all.size()); ++i)
						if (!used[i]) { g1.push_back(all[i]); b1 = b1.unionEnvelope(all[i]->getEnvelope()); used[i] = 1; ++assigned; }
					break;
				}
				if (static_cast<int>(g2.size()) + remaining == minEntries) {
					for (int i = 0; i < static_cast<int>(all.size()); ++i)
						if (!used[i]) { g2.push_back(all[i]); b2 = b2.unionEnvelope(all[i]->getEnvelope()); used[i] = 1; ++assigned; }
					break;
				}
				int choose = -1;
				double maxDiff = -std::numeric_limits<double>::max();
				for (int i = 0; i < static_cast<int>(all.size()); ++i) {
					if (used[i]) continue;
					const Envelope& eb = all[i]->getEnvelope();
					double d1 = areaEnlargement(b1, eb);
					double d2 = areaEnlargement(b2, eb);
					double diff = std::abs(d1 - d2);
					if (diff > maxDiff) { maxDiff = diff; choose = i; }
				}
				if (choose == -1) break;
				const Envelope& eb = all[choose]->getEnvelope();
				double d1 = areaEnlargement(b1, eb);
				double d2 = areaEnlargement(b2, eb);
				if (d1 < d2 || (d1 == d2 && g1.size() <= g2.size())) {
					g1.push_back(all[choose]); b1 = b1.unionEnvelope(eb);
				} else {
					g2.push_back(all[choose]); b2 = b2.unionEnvelope(eb);
				}
				used[choose] = 1; ++assigned;
			}

			// 清空inner的children
			while (inner->getChildNum() > 0) inner->popBackChildNode();
			for (auto* c : g1) inner->add(c);
			inner->setEnvelope(b1);

			newInner = new RNode(b2);
			for (auto* c : g2) newInner->add(c);
			newInner->setEnvelope(b2);
		};

		// 按 x 轴排序
		std::vector<Feature> ordered = features;
		std::sort(ordered.begin(), ordered.end(), [](const Feature& a, const Feature& b) {
			return a.getEnvelope().getMinX() < b.getEnvelope().getMinX();
		});

		// 对相同 minX 的局部段做随机乱序，以避免退化
		std::mt19937 gen(static_cast<unsigned int>(std::random_device{}()));
		size_t start = 0;
		while (start < ordered.size()) {
			double curMinX = ordered[start].getEnvelope().getMinX();
			size_t end = start + 1;
			while (end < ordered.size() && ordered[end].getEnvelope().getMinX() == curMinX)
				++end;
			if (end - start > 1) {
				std::shuffle(ordered.begin() + static_cast<std::ptrdiff_t>(start),
					ordered.begin() + static_cast<std::ptrdiff_t>(end), gen);
			}
			start = end;
		}

		// 初始化根为叶节点
		root = new RNode(ordered[0].getEnvelope());
		root->add(ordered[0]);
		// 逐个插入features（随机化后的顺序）
		for (size_t i = 1; i < ordered.size(); ++i) {
			const Feature& f = ordered[i];const Envelope& fb = f.getEnvelope();
			// 从根节点开始，选择最佳路径到达叶节点
			RNode* node = root;
			while (!node->isLeafNode()) {node = chooseSubtreeLocal(node, fb);}
			// 插入feature到叶节点
			if (static_cast<int>(node->getFeatureNum()) < maxChildren) {
				node->add(f);node->setEnvelope(node->getEnvelope().unionEnvelope(fb));updateBboxUp(node->getParent());}
			else {
				// 需要分裂
				RNode* newLeaf = nullptr;splitLeafQuadratic(node, f, newLeaf);updateBboxUp(node);
				// 向上传播分裂
				RNode* current = node;RNode* newNode = newLeaf;
				while (current->getParent() != nullptr) {
					RNode* parent = current->getParent();
					if (parent->getChildNum() < maxChildren) {
						parent->add(newNode);updateBboxUp(parent);newNode = nullptr;break;
					}
					else {
						// 分裂内部节点
						RNode* newInternal = nullptr;splitInnerQuadratic(parent, newNode, newInternal);current = parent;newNode = newInternal;
					}
				}

				// 如果分裂传播到根节点
				if (newNode != nullptr) {
					Envelope newRootBox = root->getEnvelope().unionEnvelope(newNode->getEnvelope());
					RNode* newRoot = new RNode(newRootBox);newRoot->add(root);newRoot->add(newNode);root = newRoot;}
			}
		}
		bbox = root->getEnvelope();

		return true;
	}


	void RTree::rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
		features.clear();
		if (root != nullptr)
			root->rangeQuery(rect, features);
	}
	
	bool RTree::NNQuery(double x, double y, std::vector<Feature>& features) {
		features.clear();
		// Task NNQuery
		// 使用 R 树结构做基于包围盒的最近邻过滤，只返回候选集，精炼在 hw6.cpp 中完成
		if (root == nullptr) return false;

		// 点到包围盒的最小欧氏距离
		auto minDistance2Envelope = [](double px, double py, const Envelope& e) -> double {
			double cx = px;
			if (cx < e.getMinX()) cx = e.getMinX();
			else if (cx > e.getMaxX()) cx = e.getMaxX();

			double cy = py;
			if (cy < e.getMinY()) cy = e.getMinY();
			else if (cy > e.getMaxY()) cy = e.getMaxY();

			double dx = px - cx;
			double dy = py - cy;
			return std::sqrt(dx * dx + dy * dy);
		};

		using NodeEntry = std::pair<double, RNode*>;
		struct NodeCmp {
			bool operator()(const NodeEntry& a, const NodeEntry& b) const {
				// 最小堆：距离小的优先
				return a.first > b.first;
			}
		};

		std::priority_queue<NodeEntry, std::vector<NodeEntry>, NodeCmp> pq;

		double rootMinDist = minDistance2Envelope(x, y, root->getEnvelope());
		pq.emplace(rootMinDist, root);

		// 当前找到的最小 maxDistance2Envelope 上界
		double bestMaxDist = std::numeric_limits<double>::max();

		while (!pq.empty()) {
			double nodeMinDist = pq.top().first;
			RNode* node = pq.top().second;
			pq.pop();

			// 剪枝：该节点的最小可能距离已经大于当前上界
			if (nodeMinDist > bestMaxDist) continue;

			if (node->isLeafNode()) {
				const size_t fn = node->getFeatureNum();
				for (size_t i = 0; i < fn; ++i) {
					const Feature& f = node->getFeature(i);
					double d = f.maxDistance2Envelope(x, y);
					if (d < bestMaxDist) bestMaxDist = d;
				}
			} else {
				const int cn = node->getChildNum();
				for (int i = 0; i < cn; ++i) {
					RNode* child = node->getChildNode(static_cast<size_t>(i));
					if (!child) continue;
					double childMinDist = minDistance2Envelope(x, y, child->getEnvelope());
					if (childMinDist <= bestMaxDist) {
						pq.emplace(childMinDist, child);
					}
				}
			}
		}

		// 如果没有更新过 bestMaxDist，说明树中没有要素
		if (bestMaxDist == std::numeric_limits<double>::max()) return false;

		// 使用最终上界构造查询窗口，得到候选集
		Envelope queryRect(x - bestMaxDist, x + bestMaxDist,
			y - bestMaxDist, y + bestMaxDist);
		rangeQuery(queryRect, features);
		// 注意：这里只做过滤，精确最近邻与去重在 hw6.cpp 中完成

		return !features.empty();
	}

	bool RTree::kNNQuery(double x, double y, int k, std::vector<Feature>& features) {
		// [[optional]] k最邻近几何特征查询 (k-NN)
		// 使用优先队列实现，按距离排序返回k个最近的几何特征
		features.clear();
		if (root == nullptr || k <= 0) return false;

		// 使用优先队列存储候选特征，按距离排序（最大堆，用于保持k个最小距离）
		// pair<距离, Feature>
		auto cmp = [](const std::pair<double, Feature>& a, const std::pair<double, Feature>& b) {
			return a.first < b.first; // 最大堆
		};
		std::priority_queue<std::pair<double, Feature>, 
			std::vector<std::pair<double, Feature>>, decltype(cmp)> topK(cmp);
		Point queryPoint(x, y);
		auto minDistance2Envelope = [](double px, double py, const Envelope& e) -> double {
			double cx = px;
			if (cx < e.getMinX()) cx = e.getMinX();
			else if (cx > e.getMaxX()) cx = e.getMaxX();

			double cy = py;
			if (cy < e.getMinY()) cy = e.getMinY();
			else if (cy > e.getMaxY()) cy = e.getMaxY();

			double dx = px - cx;
			double dy = py - cy;
			return std::sqrt(dx * dx + dy * dy);
		};

		// node-priority queue (min-heap by node min-distance)
		using NodeEntry = std::pair<double, RNode*>;
		struct NodeCmp { bool operator()(const NodeEntry& a, const NodeEntry& b) const { return a.first > b.first; } };
		std::priority_queue<NodeEntry, std::vector<NodeEntry>, NodeCmp> nodePQ;
		double rootMin = minDistance2Envelope(x, y, root->getEnvelope());
		nodePQ.emplace(rootMin, root);std::unordered_set<const Geometry*> seen;seen.reserve(256);
		double bestKDist = std::numeric_limits<double>::infinity(); // current worst distance in top-k
		while (!nodePQ.empty()) {
			auto entry = nodePQ.top(); nodePQ.pop();double nodeMinDist = entry.first;RNode* node = entry.second;
			// pruning: if minimal possible distance to this node is greater than current k-th best
			if (nodeMinDist > bestKDist) break;
			if (node->isLeafNode()) {
				const size_t fn = node->getFeatureNum();
				for (size_t i = 0; i < fn; ++i) {
					const Feature& f = node->getFeature(i);const Geometry* geom = f.getGeom();if (!geom) continue;
					if (!seen.insert(geom).second) continue; // 去重
					double d = geom->distance(&queryPoint);
					if (static_cast<int>(topK.size()) < k) {topK.push(std::make_pair(d, f));
						if (static_cast<int>(topK.size()) == k) bestKDist = topK.top().first;}
					else if (d < topK.top().first) {
						topK.pop();topK.push(std::make_pair(d, f));bestKDist = topK.top().first;}}}
			else {const int cn = node->getChildNum();
				for (int i = 0; i < cn; ++i) {
					RNode* child = node->getChildNode(static_cast<size_t>(i));
					if (!child) continue;
					double childMin = minDistance2Envelope(x, y, child->getEnvelope());
					// only consider child if it may contain closer objects than current worst
					if (childMin <= bestKDist) nodePQ.emplace(childMin, child);}}
		}
		if (topK.empty()) return false;
		// 从优先队列中提取结果（需要反转顺序）
		std::vector<std::pair<double, Feature>> temp;while (!topK.empty()) { temp.push_back(topK.top()); topK.pop(); }
		// 按距离从小到大排序并输出
		for (auto it = temp.rbegin(); it != temp.rend(); ++it) {features.push_back(it->second);}
		return !features.empty();}

	void RTree::spatialJoin(const std::vector<Feature>& A,
		const std::vector<Feature>& B,double dist,std::vector<std::pair<Feature, Feature>>& result) {
		/*
		 * Spatial Join 基于距离的空间关联
		 * Filter 步骤：对集合 A 中每个几何对象 a，扩展其包围盒 dist 距离，
		 *   使用 R 树在集合 B 中做范围查询，得到候选集合；
		 * Refine 步骤：对候选集合计算精确几何距离，保留距离 <= dist 的配对；
		 * 正确性：任意满足距离约束的配对，其几何对象的包围盒必与扩展后的
		 *   包围盒相交，因此不会被 Filter 步骤丢弃。
		 */
		result.clear();
		if (A.empty() || B.empty())
			return;
		// 针对 B 构建局部 R 树索引
		RTree treeB(this->maxChildren);
		if (!treeB.constructTree(B)) {
			return;
		}
		std::vector<Feature> candidates;candidates.reserve(64);
		for (const Feature& a : A) {
			candidates.clear(); // 每个 a 单独收集候选，避免累加
			const Envelope& eb = a.getEnvelope();
			Envelope searchEnv(eb.getMinX() - dist, eb.getMaxX() + dist,
			eb.getMinY() - dist, eb.getMaxY() + dist);
			// Filter：用 treeB 在 B 中做范围查询
			treeB.rangeQuery(searchEnv, candidates);
			// Refine：使用精确距离，并基于 Geometry* 去重
			std::unordered_set<const Geometry*> seenB;
			seenB.reserve(candidates.size() * 2);
			for (const Feature& b : candidates) {
				const Geometry* gb = b.getGeom();const Geometry* ga = a.getGeom();
				if (!ga || !gb) continue;
				if (!seenB.insert(gb).second) continue; // 同一几何只检查一次
				double d = ga->distance(gb);
				if (d <= dist) {result.emplace_back(a, b);}
			}
		}
	}

} // namespace hw6

