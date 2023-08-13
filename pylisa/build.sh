make clean

echo "Make linux shared library"
export OS_OPTION="linux"
make

echo "Test shared library"
python test_mie_wrapper.py

echo "Make windows dll"
export OS_OPTION="windows"
make

echo "Make for darwin"
export OS_OPTION="darwin"
make