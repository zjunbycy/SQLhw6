#include "RTree.h"

namespace hw6 {

	// RNode 实现
	void RNode::add(RNode* child) {
		children[childrenNum] = child;
		child->parent = this;
		++childrenNum;
	}

	void RNode::remove(const Feature& f) {
		auto where = [&]() {
			for (auto itr = features.begin(); itr != features.end(); ++itr)
				if (itr->getName() == f.getName())
					return itr;
			}();
			features.erase(where);
			if (features.empty())
				features.shrink_to_fit(); // free memory unused but allocated
	}

	void RNode::remove(RNode* child) {
		for (int i = 0; i < childrenNum; ++i)
			if (children[i] == child) {
				--childrenNum;
				std::swap(children[i], children[childrenNum]);
				children[childrenNum] = nullptr;
				break;
			}
	}

	Feature RNode::popBackFeature() {
		auto ret = features.back();
		features.pop_back();
		return ret;
	}

	RNode* RNode::popBackChildNode() {
		--childrenNum;
		auto ret = children[childrenNum];
		children[childrenNum] = nullptr;
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
			for (int i = 0; i < childrenNum; ++i) {
				children[i]->rangeQuery(rect, features);
			}
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
		// 按x轴（包围盒minX）顺序插入，节点选择面积增量最小原则，溢出使用二次分裂，种子为最左/最右
		if (features.empty())
			return false;

		// 辅助：面积增量
		auto areaEnlargement = [](const Envelope& nodeBox, const Envelope& box) {
			Envelope u = nodeBox.unionEnvelope(box);
			return u.getArea() - nodeBox.getArea();
		};

		// 辅助：从当前内部节点选择子树（面积增量最小，若相同选择面积更小）
		auto chooseSubtreeLocal = [&](RNode* node, const Envelope& box) -> RNode* {
			RNode* best = nullptr;
			double minEnl = std::numeric_limits<double>::max();
			double minArea = std::numeric_limits<double>::max();
			for (int i = 0; i < node->getChildNum(); ++i) {
				RNode* c = node->getChildNode(static_cast<size_t>(i));
				const Envelope& cb = c->getEnvelope();
				double enl = areaEnlargement(cb, box);
				double area = cb.getArea();
				if (enl < minEnl || (enl == minEnl && area < minArea)) {
					minEnl = enl;
					minArea = area;
					best = c;
				}
			}
			return best;
		};

		// 辅助：向上更新bbox
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
			g1.push_back(all[leftIdx]);
			g2.push_back(all[rightIdx]);
			Envelope b1 = all[leftIdx].getEnvelope();
			Envelope b2 = all[rightIdx].getEnvelope();

			std::vector<char> used(all.size(), 0);
			used[leftIdx] = used[rightIdx] = 1;

			int minEntries = std::max(1, maxChildren / 2);
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
					double d1 = areaEnlargement(b1, eb);
					double d2 = areaEnlargement(b2, eb);
					double diff = std::abs(d1 - d2);
					if (diff > maxDiff) { maxDiff = diff; choose = i; }
				}
				if (choose == -1) break;
				const Envelope& eb = all[choose].getEnvelope();
				double d1 = areaEnlargement(b1, eb);
				double d2 = areaEnlargement(b2, eb);
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
			g1.push_back(all[leftIdx]);
			g2.push_back(all[rightIdx]);
			Envelope b1 = all[leftIdx]->getEnvelope();
			Envelope b2 = all[rightIdx]->getEnvelope();

			std::vector<char> used(all.size(), 0);
			used[leftIdx] = used[rightIdx] = 1;

			int minEntries = std::max(1, maxChildren / 2);
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

		// 按x轴顺序插入（minX升序）
		std::vector<Feature> ordered = features;
		std::sort(ordered.begin(), ordered.end(), [](const Feature& a, const Feature& b) {
			return a.getEnvelope().getMinX() < b.getEnvelope().getMinX();
		});

		// 初始化根为叶节点
		root = new RNode(ordered[0].getEnvelope());
		root->add(ordered[0]);

		// 逐个插入features
		for (size_t i = 1; i < ordered.size(); ++i) {
			const Feature& f = ordered[i];
			const Envelope& fb = f.getEnvelope();

			// 从根节点开始，选择最佳路径到达叶节点
			RNode* node = root;
			while (!node->isLeafNode()) {
				node = chooseSubtreeLocal(node, fb);
			}

			// 插入feature到叶节点
			if (static_cast<int>(node->getFeatureNum()) < maxChildren) {
				node->add(f);
				node->setEnvelope(node->getEnvelope().unionEnvelope(fb));
				updateBboxUp(node->getParent());
			}
			else {
				// 需要分裂
				RNode* newLeaf = nullptr;
				splitLeafQuadratic(node, f, newLeaf);
				updateBboxUp(node);

				// 向上传播分裂
				RNode* current = node;
				RNode* newNode = newLeaf;
				while (current->getParent() != nullptr) {
					RNode* parent = current->getParent();
					if (parent->getChildNum() < maxChildren) {
						parent->add(newNode);
						updateBboxUp(parent);
						newNode = nullptr;
						break;
					}
					else {
						// 分裂内部节点
						RNode* newInternal = nullptr;
						splitInnerQuadratic(parent, newNode, newInternal);
						current = parent;
						newNode = newInternal;
					}
				}

				// 如果分裂传播到根节点
				if (newNode != nullptr) {
					Envelope newRootBox = root->getEnvelope().unionEnvelope(newNode->getEnvelope());
					RNode* newRoot = new RNode(newRootBox);
					newRoot->add(root);
					newRoot->add(newNode);
					root = newRoot;
				}
			}
		}

		// 更新整棵树的包围盒
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
		/* TODO */
		if(root == nullptr){return false;}
		// filter step
		// (使用maxDistance2Envelope函数，获得查询点到几何对象包围盒的最短的最大距离，然后区域查询获得候选集)
		// 首先找到包含查询点的叶节点，或最近的叶节点
		RNode* leafNode = root->pointInLeafNode(x, y);
		double minMaxDist = std::numeric_limits<double>::max();
		if (leafNode != nullptr && leafNode->getFeatureNum() > 0) {
			// 在叶节点中找到最短的最大距离
			for (size_t i = 0; i < leafNode->getFeatureNum(); ++i) {
				double dist = leafNode->getFeature(i).maxDistance2Envelope(x, y);
				minMaxDist = std::min(minMaxDist, dist);
			}
		}
		else {
			// 查询点不在任何叶节点内，遍历所有叶节点找最小maxDistance
			std::queue<RNode*> queue;
			queue.push(root);
			while (!queue.empty()) {
				RNode* node = queue.front();
				queue.pop();
				if (node->isLeafNode()) {
					for (size_t i = 0; i < node->getFeatureNum(); ++i) {
						double dist = node->getFeature(i).maxDistance2Envelope(x, y);
						minMaxDist = std::min(minMaxDist, dist);
					}
				}
				else {
					for (int i = 0; i < node->getChildNum(); ++i) {
						queue.push(node->getChildNode(i));
					}
				}
			}
		}
		Envelope queryRect(x - minMaxDist, x + minMaxDist, y - minMaxDist, y + minMaxDist);
		rangeQuery(queryRect, features);
		// 注意R树邻近查询仅返回候选集，精炼步在hw6的NNQuery中完成

		return !features.empty();
	}


} // namespace hw6

