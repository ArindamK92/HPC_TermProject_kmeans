g++ -o3 -o op_kmean_base main.cpp
g++ -o3 -o op_kmean_reg register_test.cpp
g++ -o3 -o op_kmean_regBlock register_test_blockwise.cpp

./op_kmean_base <n> <m> <k> <itr_max>
./op_kmean_reg <n> <m> <k> <itr_max>
./op_kmean_regBlock <n> <m> <k> <itr_max>


./op_kmean_base 33554432 32 16 5