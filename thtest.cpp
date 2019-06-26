#include <iostream>
#include <thread>
#include <future>
#include <functional>
#include <chrono>
#include <vector>

using namespace std;

void doPrintVals(promise<vector<int> > &p) {
	vector<int> ls;
	thread* t1;
	auto lambda = ([&ls]() void {
		for(int i=0;i<5;i++) {
			ls.push_back(i);
			cout<<"child: "<<i<<endl;
	   		std::this_thread::sleep_for (std::chrono::milliseconds(300));
		}
	});
	t1 = new thread(lambda);
	t1.join();
	cout<<"ls size: "<<ls.size()<<endl;
	p.set_value(ls);
}

int main() {

	promise<vector<int> > p;
	thread* t2;
	t2 = new thread(doPrintVals, ref(p));
	for(int i=0;i<10;i++) {
		cout<<"main: "<<i<<endl;
   		std::this_thread::sleep_for (std::chrono::milliseconds(100));
	}

	t2->join();
	future<vector<int> > f = p.get_future();
	vector<int> result = f.get();
	cout<<"result_size: "<<result.size()<<endl;
}
