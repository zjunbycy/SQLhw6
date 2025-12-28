#include "Common.h"
#include "RTree.h"
#include <ctime>
#include <set>
#include <utility>
#include "CMakeIn.h"

using namespace hw6;

extern int mode;
extern std::vector<Geometry*> readGeom(const char* filename);
extern std::vector<std::string> readName(const char* filename);
extern void transformValue(double& res, const char* format);
extern void wrongMessage(Envelope e1, Envelope e2, bool cal);
extern void wrongMessage(const Point& pt1, const Point& pt2, double dis,
    double res);
extern void wrongMessage(Envelope e1, Envelope e2, Envelope cal, Envelope res);

// 声明外部变量用于Spatial Join测试
extern std::vector<Feature> features;
extern std::vector<Feature> roads;

namespace hw6 {

    void RTree::test(int t) {
        using namespace std;

        std::cout << "*********************Start*********************" << std::endl;
        if (t == TEST1) {
            std::cout << "TEST1: Envelope Contain, Intersect, and Union" << endl;

            int failedCase = 0;
            Envelope e1(-1, 1, -1, 1);
            vector<Envelope> tests;
            tests.push_back(Envelope(-0.5, 0.5, -0.5, 0.5));
            tests.push_back(Envelope(-0.5, 0.5, 0.5, 1.5));
            tests.push_back(Envelope(0.5, 1.5, -0.5, 0.5));
            tests.push_back(Envelope(-1.5, -0.5, -1.5, -0.5));
            tests.push_back(Envelope(-2, -1, -0.5, 0.5));
            tests.push_back(Envelope(1, 1.5, 1, 1.5));
            tests.push_back(Envelope(-2, -1.5, -0.5, 0.5));
            tests.push_back(Envelope(-0.5, 0.5, 1.5, 2));
            tests.push_back(Envelope(-2, -1.5, 0.5, 1.5));
            tests.push_back(Envelope(0.5, 1.5, 1.5, 2));

            for (size_t i = 0; i < tests.size(); ++i) {
                if (e1.contain(tests[i]) != (i == 0)) {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], (i != 0));
                }
                if (tests[i].contain(e1) == true) {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, true);
                }
            }
            cout << "Envelope Contain: " << tests.size() * 2 - failedCase << " / "
                << tests.size() * 2 << " tests are passed" << endl;

            failedCase = 0;
            for (size_t i = 0; i < tests.size(); ++i) {
                if (e1.intersect(tests[i]) != (i < 6)) {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], (i < 6));
                }
                if (tests[i].intersect(e1) != (i < 6)) {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, (i < 6));
                }
            }
            cout << "Envelope Intersect: " << tests.size() * 2 - failedCase << " / "
                << tests.size() * 2 << " tests are passed" << endl;

            failedCase = 0;
            vector<Envelope> results;
            results.push_back(Envelope(-1, 1, -1, 1));
            results.push_back(Envelope(-1, 1, -1, 1.5));
            results.push_back(Envelope(-1, 1.5, -1, 1));
            results.push_back(Envelope(-1.5, 1, -1.5, 1));
            results.push_back(Envelope(-2, 1, -1, 1));
            results.push_back(Envelope(-1, 1.5, -1, 1.5));
            results.push_back(Envelope(-2, 1, -1, 1));
            results.push_back(Envelope(-1, 1, -1, 2));
            results.push_back(Envelope(-2, 1, -1, 1.5));
            results.push_back(Envelope(-1, 1.5, -1, 2));
            for (size_t i = 0; i < tests.size(); ++i) {
                if (e1.unionEnvelope(tests[i]) != results[i]) {
                    failedCase += 1;
                    wrongMessage(e1, tests[i], e1.unionEnvelope(tests[i]),
                        results[i]);
                }
                if (tests[i].unionEnvelope(e1) != results[i]) {
                    failedCase += 1;
                    wrongMessage(tests[i], e1, e1.unionEnvelope(tests[i]),
                        results[i]);
                }
            }
            cout << "Envelope Union: " << tests.size() * 2 - failedCase << " / "
                << tests.size() * 2 << " tests are passed" << endl;
        }
        else if (t == TEST2) {
            cout << "TEST2: Distance between Point and LineString" << endl;

            vector<Point> points;
            points.push_back(Point(0, 0));
            points.push_back(Point(10, 10));
           	LineString line(points);

            points.push_back(Point(-10, -10));
            points.push_back(Point(20, 20));
            points.push_back(Point(5, 5));
            points.push_back(Point(10, 0));
            points.push_back(Point(10, -10));
            points.push_back(Point(0, 10));
            points.push_back(Point(0, 20));
            points.push_back(Point(20, 0));

            double dists[] = { 0,       0,       14.1421, 14.1421, 0,
							  7.07107, 14.1421, 7.07107, 14.1421, 14.1421 };

            int failedCase = 0;
            for (size_t i = 0; i < points.size(); ++i) {
                double dist = points[i].distance(&line);
                if (fabs(dist - dists[i]) > 0.0001) {
                    failedCase += 1;
                    cout << "Your answer is " << dist << " for test between ";
                    line.print();
                    cout << " and ";
                    points[i].print();
                    cout << ", but the answer is " << dists[i] << endl;
                }
            }
            cout << "Distance between Point and LineString: "
                << points.size() - failedCase << " / " << points.size()
                << " tests are passed" << endl;
        }
        else if (t == TEST3) {
            cout << "TEST3: Distance between Point and Polygon" << endl;

            vector<Point> points;
            points.push_back(Point(5, 0));
            points.push_back(Point(3, 6));
            points.push_back(Point(2, 4));
            points.push_back(Point(-2, 4));
            points.push_back(Point(-3, 5));
            points.push_back(Point(-5, 0));
            points.push_back(Point(0, -3));
            points.push_back(Point(5, 0));
            LineString line(points);
            Polygon poly(line);

            points.clear();
            points.push_back(Point(5, 4));
            points.push_back(Point(3, 4));
            points.push_back(Point(0, 4));
            points.push_back(Point(-3, 4));
            points.push_back(Point(-5, 4));
            points.push_back(Point(5, 5));
            points.push_back(Point(3, 5));
            points.push_back(Point(0, 5));
            points.push_back(Point(-3, 5));
            points.push_back(Point(0, 0));

            double dists[] = { 1.26491, 0, 0, 0, 1.48556, 1.58114, 0, 1, 0, 0 };

            int failedCase = 0;
            for (size_t i = 0; i < points.size(); ++i) {
                double dist = points[i].distance(&poly);
                if (fabs(dist - dists[i]) > 0.00001) {
                    failedCase += 1;
                    cout << "Your answer is " << dist << " for test between ";
                    poly.print();
                    cout << " and ";
                    points[i].print();
                    cout << ", but the answer is " << dists[i] << endl;
                }
            }
            cout << "Distance between Point and Polygon: "
                << points.size() - failedCase << " / " << points.size()
                << " tests are passed" << endl;
        }
        else if (t == TEST4) {
            cout << "TEST4: RTree Construction" << endl;
            int ncase, cct;
            ncase = cct = 2;

            RTree rtree(8);
            vector<Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/station");
            vector<Feature> features;

            for (size_t i = 0; i < geom.size(); ++i)
                features.push_back(Feature("", geom[i]));

            rtree.constructTree(features);

            int height, interiorNum, leafNum;
            rtree.countHeight(height);
            rtree.countNode(interiorNum, leafNum);

            if (!(height == 4 && interiorNum == 18 && leafNum == 86)) {
                cout << "Case 1: "
                    << "Your answer is height: " << height
                    << ", interiorNum: " << interiorNum << ", leafNum: " << leafNum
                    << ". One possible answer is height: 4, interiorNum: "
                    "18, "
                    "leafNum: 86\n";
                --cct;
            }

            features.clear();
            for (size_t i = 0; i < geom.size(); ++i)
                delete geom[i];
            geom.clear();

            vector<Geometry*> geom2 = readGeom(PROJ_SRC_DIR "/data/highway");
            vector<Feature> features2;
            RTree rtree2(8);

            for (size_t i = 0; i < geom2.size(); ++i)
                features2.push_back(Feature("", geom2[i]));

            rtree2.constructTree(features2);

            int height2, interiorNum2, leafNum2;
            rtree2.countHeight(height2);
            rtree2.countNode(interiorNum2, leafNum2);

            if (!(height2 == 6 && interiorNum2 == 484 && leafNum2 == 2305)) {
                cout << "Case 2: "
                    << "Your answer is height: " << height2
                    << ", interiorNum: " << interiorNum2
                    << ", leafNum: " << leafNum2
                    << ". One possible answer is height: 6, interiorNum: "
                    "484, leafNum: 2305\n";
                --cct;
            }

            features2.clear();
            for (size_t i = 0; i < geom2.size(); ++i)
                delete geom2[i];
            geom2.clear();

            //cout << "RTree Construction: " << cct << " / " << ncase
            //     << " tests are passed" << endl;
        }
        else if (t == TEST5) {
            // TEST5: Spatial Join (Station-Road) + correctness/perf similar to QuadTreeTest
            cout << "TEST5: Spatial Join with Correctness Verification" << endl;
            cout << "=================================================" << endl;

            // 读取 station 和 highway 数据
            vector<Geometry*> geomA = readGeom(PROJ_SRC_DIR "/data/station");
            vector<Geometry*> geomB = readGeom(PROJ_SRC_DIR "/data/highway");

            vector<Feature> A, B;
            A.reserve(geomA.size());
            B.reserve(geomB.size());

            for (size_t i = 0; i < geomA.size(); ++i)
                A.push_back(Feature("", geomA[i]));
            for (size_t i = 0; i < geomB.size(); ++i)
                B.push_back(Feature("", geomB[i]));

            cout << "Dataset: " << A.size() << " stations x " << B.size() << " roads" << endl;

            // 执行 Spatial Join
            RTree rtree(8);
            vector<pair<Feature, Feature>> result;
            double joinDist = 0.001;

            clock_t start = clock();
            rtree.spatialJoin(A, B, joinDist, result);
            clock_t end = clock();

            double rtree_time = (end - start) / 1000.0;

            cout << "Distance threshold: " << joinDist << endl;
            cout << "Pairs found: " << result.size() << endl;
            cout << "RTree time: " << rtree_time << "s" << endl;

            // === 性能对比：暴力方法 vs RTree ===
            cout << "\n--- Performance Comparison ---" << endl;

            int testA = (100 < (int)A.size()) ? 100 : (int)A.size();
            int testB = (100 < (int)B.size()) ? 100 : (int)B.size();

            clock_t brute_start = clock();
            int brute_sample_result = 0;
            for (int i = 0; i < testA; ++i) {
                for (int j = 0; j < testB; ++j) {
                    if (A[i].getGeom()->distance(B[j].getGeom()) <= joinDist) {
                        brute_sample_result++;
                    }
                }
            }
            clock_t brute_end = clock();
            double brute_sample_time = (brute_end - brute_start) / 1000.0;

            double total_comparisons = (double)A.size() * B.size();
            double sample_comparisons = (double)testA * testB;
            double estimated_brute_time = brute_sample_time * (total_comparisons / sample_comparisons);

            cout << "Brute force sample (" << testA << "x" << testB << "): " 
                 << brute_sample_time << "s (" << brute_sample_result << " pairs)" << endl;
            cout << "Estimated full brute force time: " << estimated_brute_time << "s" << endl;
            cout << "Speedup: " << estimated_brute_time / rtree_time << "x" << endl;

            // === 正确性验证 ===
            cout << "\n--- Correctness Verification ---" << endl;

            int sampleSize = (50 < (int)result.size()) ? 50 : (int)result.size();
            int validCount = 0;
            for (int i = 0; i < sampleSize; ++i) {
                double dist = result[i].first.getGeom()->distance(result[i].second.getGeom());
                if (dist <= joinDist) {
                    validCount++;
                }
            }

            cout << "Distance constraint check: " << validCount << "/" << sampleSize 
                 << (validCount == sampleSize ? " PASS" : " FAIL") << endl;

            std::set<std::pair<size_t, size_t>> uniquePairs;
            for (const auto& p : result) {
                size_t addr1 = (size_t)p.first.getGeom();
                size_t addr2 = (size_t)p.second.getGeom();
                uniquePairs.insert(std::make_pair(addr1, addr2));
            }

            cout << "Deduplication check: " << uniquePairs.size() << "/" << result.size()
                 << (uniquePairs.size() == result.size() ? " PASS" : " FAIL") << endl;

            int rtreeCount = 0;
            for (const auto& p : result) {
                bool aFound = false, bFound = false;
                for (int i = 0; i < testA && !aFound; ++i) {
                    if (A[i].getGeom() == p.first.getGeom()) aFound = true;
                }
                for (int j = 0; j < testB && !bFound; ++j) {
                    if (B[j].getGeom() == p.second.getGeom()) bFound = true;
                }
                if (aFound && bFound) rtreeCount++;
            }

            cout << "Brute force validation (" << testA << "x" << testB << "): "
                 << "brute=" << brute_sample_result << ", rtree=" << rtreeCount
                 << (brute_sample_result == rtreeCount ? " PASS" : " FAIL") << endl;

            cout << "\nCorrectness verification: "
                 << ((validCount == sampleSize) && 
                     (uniquePairs.size() == result.size()) && 
                     (brute_sample_result == rtreeCount) ? "ALL PASS" : "SOME FAIL")
                 << endl;
            cout << "=================================================\n" << endl;

            // 清理几何内存
            for (size_t i = 0; i < geomA.size(); ++i)
                delete geomA[i];
            for (size_t i = 0; i < geomB.size(); ++i)
                delete geomB[i];
            geomA.clear();
            geomB.clear();
        }
        else if (t == TEST6) {
            // TEST6: k-NN 查询测试
            cout << "TEST6: k-NN Query Test" << endl;
			
			RTree rtree(8);
			vector<Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/station");
			vector<string> name = readName(PROJ_SRC_DIR "/data/station");
			vector<Feature> features;
			
			for (size_t i = 0; i < geom.size(); ++i)
				features.push_back(Feature(name[i], geom[i]));
			
			rtree.constructTree(features);
			
			// 测试点（使用数据集中心附近的点）
			Envelope bbox = rtree.getEnvelope();
			double cx = (bbox.getMinX() + bbox.getMaxX()) / 2;
			double cy = (bbox.getMinY() + bbox.getMaxY()) / 2;
			
			cout << "Query point: (" << cx << ", " << cy << ")" << endl;
			
			// 测试不同的k值
			int kValues[] = { 1, 3, 5, 10 };
			for (int k : kValues) {
				vector<Feature> knnResult;
				rtree.kNNQuery(cx, cy, k, knnResult);
				
				cout << "k=" << k << ", found " << knnResult.size() << " neighbors:" << endl;
				Point queryPt(cx, cy);
				for (size_t i = 0; i < knnResult.size(); ++i) {
					double dist = knnResult[i].getGeom()->distance(&queryPt);
					cout << "  " << (i+1) << ". " << knnResult[i].getName() 
						 << ", distance: " << dist << endl;
				}
			}
			
			// 清理内存
			for (auto* g : geom) delete g;
        }
        else if (t == TEST8) {
            cout << "TEST8: RTreeAnalysis" << endl;
            analyse();
        }

        cout << "**********************End**********************" << endl;
    }

    void forConstCapAnalyseRTree(const std::vector<Feature>& features, int childNum, int maxNum, int step) {
        if (childNum <= maxNum) {
            RTree rtree(childNum);
            // TODO 与四叉树进行比较

            // 构造R树，输出R树的节点数目和高度
            clock_t start_time = clock();
            rtree.constructTree(features);
            clock_t end_time = clock();

            int height = 0, interiorNum = 0, leafNum = 0;
            rtree.countHeight(height);
            rtree.countNode(interiorNum, leafNum);

            std::cout << "MaxChildren " << childNum << "\n";
            std::cout << "Height: " << height
                << " \tInterior node number: " << interiorNum
                << " \tLeaf node number: " << leafNum << "\n";
            std::cout << "Construction time: "
                << (end_time - start_time) / 1000.0 << "s" << std::endl;

            // 评估100000次随机最邻近几何特征查询的时间
            double x, y;
            std::vector<Feature> candidateFeatures;
            start_time = clock();
            for (int i = 0; i < 100000; ++i) {
                x = -((rand() % 225) / 10000.0 + 73.9812);
                y = (rand() % 239) / 10000.0 + 40.7247;
                rtree.NNQuery(x, y, candidateFeatures);
                // refine step 在 hw6.cpp 中完成，这里只统计索引过滤部分
                candidateFeatures.clear();
            }
            end_time = clock();
            std::cout << "NNQuery time: "
                << (end_time - start_time) / 1000.0 << "s" << std::endl
                << std::endl;

            forConstCapAnalyseRTree(features, childNum + step, maxNum, step);
        }
    }

    void RTree::analyse() {
        using namespace std;

        vector<Feature> features;
        vector<Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/taxi");
        vector<string> name = readName(PROJ_SRC_DIR "/data/taxi");

        features.clear();
        features.reserve(geom.size());
        for (size_t i = 0; i < geom.size(); ++i)
            features.push_back(Feature(name[i], geom[i]));

        cout << "taxi number: " << geom.size() << endl;

        srand(time(nullptr));

        /*TODO:实现forConstCapAnalyseRTree */
        // 与 QuadTree::analyse 类似，这里在 [70, 200] 范围内
        // 变化 RTree 的 maxChildren 参数，分析其高度、节点数和 NNQuery 时间
        forConstCapAnalyseRTree(features, 70, 200, 10);
    }


} // namespace hw6
