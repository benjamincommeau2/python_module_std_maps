#g++ -shared -fPIC -I/usr/include/python3.6m -I ~/boost/boost_1_72_0 test.cpp -o test.o -fmax-errors=$1
g++ -I/usr/include/python3.6m test.cpp -o test.o -fmax-errors=$1
#g++ -I/usr/include/python2.7 -lboost_python -lboost_system test.cpp -o test.o -fmax-errors=$1
#g++ test.cpp -o test.o -fmax-errors=$1
