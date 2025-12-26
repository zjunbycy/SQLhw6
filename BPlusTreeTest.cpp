#include "BPlusTree.h"
#include "Common.h"
#include <ctime>
#include "CMakeIn.h"

using namespace hw6;

extern int mode;
extern std::vector<Geometry*> readGeom(const char* filename);
extern std::vector<std::string> readName(const char* filename);

namespace hw6 {

	void BPlusTree::test(int t) {
		using namespace std;

		cout << "*********************Start BPlusTree Test*********************" << endl;

		// TEST1: Envelope 测试（通用几何测试）
		if (t == TEST1) {
			cout << "TEST1: Envelope Contain, Intersect, and Union" << endl;
			cout << "(This is a geometry test, not specific to B+Tree)" << endl;
			
			// 这些测试与树无关，可以直接运行
			int failedCase = 0;
			Envelope e1(-1, 1, -1, 1);
			vector<Envelope> tests;
			tests.push_back(Envelope(-0.5, 0.5, -0.5, 0.5));
			tests.push_back(Envelope(-0.5, 0.5, 0.5, 1.5));
			tests.push_back(Envelope(0.5, 1.5, -0.5, 0.5));
			tests.push_back(Envelope(-1.5, -0.5, -1.5, -0.5));
			tests.push_back(Envelope(-2, -1, -0.5, 0.5));
			tests.push_back(Envelope(1, 1.5, 1, 1.5));

			// Contain 测试
			for (size_t i = 0; i < tests.size(); ++i) {
				if (e1.contain(tests[i]) != (i == 0)) failedCase++;
				if (tests[i].contain(e1) == true) failedCase++;
			}
			cout << "Envelope Contain: " << tests.size() * 2 - failedCase 
				<< " / " << tests.size() * 2 << " passed" << endl;

			// Intersect 测试
			failedCase = 0;
			for (size_t i = 0; i < tests.size(); ++i) {
				if (e1.intersect(tests[i]) != (i < 6)) failedCase++;
				if (tests[i].intersect(e1) != (i < 6)) failedCase++;
			}
			cout << "Envelope Intersect: " << tests.size() * 2 - failedCase 
				<< " / " << tests.size() * 2 << " passed" << endl;
		}
		// TEST2: Point-LineString 距离测试
		else if (t == TEST2) {
			cout << "TEST2: Distance between Point and LineString" << endl;
			cout << "(This is a geometry test, not specific to B+Tree)" << endl;

			vector<Point> points;
			points.push_back(Point(0, 0));
			points.push_back(Point(10, 10));
			LineString line(points);

			points.push_back(Point(-10, -10));
			points.push_back(Point(5, 5));
			
			int passed = 0;
			if (fabs(points[0].distance(&line) - 0.0) < 0.001) passed++;
			if (fabs(points[1].distance(&line) - 0.0) < 0.001) passed++;
			
			cout << "Point-LineString distance: " << passed << " / 2 tests passed" << endl;
		}
		// TEST3: Point-Polygon 距离测试
		else if (t == TEST3) {
			cout << "TEST3: Distance between Point and Polygon" << endl;
			cout << "(This is a geometry test, not specific to B+Tree)" << endl;

			vector<Point> points;
			points.push_back(Point(5, 0));
			points.push_back(Point(0, 5));
			points.push_back(Point(-5, 0));
			points.push_back(Point(0, -5));
			points.push_back(Point(5, 0));
			LineString line(points);
			Polygon poly(line);

			Point testPt(0, 0);
			double dist = testPt.distance(&poly);
			cout << "Distance from (0,0) to polygon: " << dist << endl;
			cout << "(Expected: 0.0 if point is inside)" << endl;
		}
		// TEST4: B+树构建测试
		else if (t == TEST4) {
			cout << "TEST4: BPlusTree Construction" << endl;

			BPlusTree bptree(20, 1024);
			vector<Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/station");
			vector<Feature> features;

			for (size_t i = 0; i < geom.size(); ++i)
				features.push_back(Feature("", geom[i]));

			clock_t start = clock();
			bptree.constructTree(features);
			clock_t end = clock();

			int height, interiorNum, leafNum;
			bptree.countHeight(height);
			bptree.countNode(interiorNum, leafNum);

			cout << "Construction time: " << (end - start) / 1000.0 << "s" << endl;
			cout << "Height: " << height
				<< " \tInterior nodes: " << interiorNum
				<< " \tLeaf nodes: " << leafNum << endl;

			// 清理
			for (size_t i = 0; i < geom.size(); ++i)
				delete geom[i];
		}
		// TEST5: 空间连接测试
		else if (t == TEST5) {
			cout << "TEST5: BPlusTree Spatial Join" << endl;

			vector<Geometry*> geomA = readGeom(PROJ_SRC_DIR "/data/station");
			vector<Geometry*> geomB = readGeom(PROJ_SRC_DIR "/data/highway");

			vector<Feature> A, B;
			for (size_t i = 0; i < geomA.size(); ++i)
				A.push_back(Feature("", geomA[i]));
			for (size_t i = 0; i < geomB.size(); ++i)
				B.push_back(Feature("", geomB[i]));

			BPlusTree bptree(20, 1024);
			vector<pair<Feature, Feature>> result;
			double joinDist = 0.001;

			clock_t start = clock();
			bptree.spatialJoin(A, B, joinDist, result);
			clock_t end = clock();

			cout << "Spatial Join distance = " << joinDist
				<< ", pairs found = " << result.size() << endl;
			cout << "Join time: " << (end - start) / 1000.0 << "s" << endl;

			// 清理
			for (size_t i = 0; i < geomA.size(); ++i) delete geomA[i];
			for (size_t i = 0; i < geomB.size(); ++i) delete geomB[i];
		}
		// TEST8: 性能分析
		else if (t == TEST8) {
			cout << "TEST8: BPlusTree Performance Analysis" << endl;
			analyse();
		}
		// 其他测试
		else {
			cout << "TEST" << t << ": Not implemented for B+Tree" << endl;
			cout << "Available tests: 1, 2, 3, 4, 5, 8" << endl;
		}

		cout << "**********************End BPlusTree Test**********************" << endl;
	}

	void BPlusTree::analyse() {
		using namespace std;

		vector<Feature> features;
		vector<Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/taxi");
		vector<string> name = readName(PROJ_SRC_DIR "/data/taxi");

		features.reserve(geom.size());
		for (size_t i = 0; i < geom.size(); ++i)
			features.push_back(Feature(name[i], geom[i]));

		cout << "=== BPlusTree Analysis ===" << endl;
		cout << "Taxi number: " << geom.size() << endl;

		srand(time(nullptr));

		for (int cap = 70; cap <= 200; cap += 10) {
			BPlusTree bptree(cap, 1024);

			clock_t start = clock();
			bptree.constructTree(features);
			clock_t end = clock();

			int height, interiorNum, leafNum;
			bptree.countHeight(height);
			bptree.countNode(interiorNum, leafNum);

			cout << "\nCapacity: " << cap << endl;
			cout << "Height: " << height
				<< " \tInterior: " << interiorNum
				<< " \tLeaf: " << leafNum << endl;
			cout << "Construction time: " << (end - start) / 1000.0 << "s" << endl;

			// 测试NN查询性能
			vector<Feature> candidateFeatures;
			start = clock();
			for (int i = 0; i < 100000; ++i) {
				double x = -((rand() % 225) / 10000.0 + 73.9812);
				double y = (rand() % 239) / 10000.0 + 40.7247;
				bptree.NNQuery(x, y, candidateFeatures);
				candidateFeatures.clear();
			}
			end = clock();
			cout << "NNQuery time: " << (end - start) / 1000.0 << "s" << endl;
		}
	}

} // namespace hw6