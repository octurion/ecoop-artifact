ROOT_PATH=op2orig/apps/c/build

COUNT=20

cd $ROOT_PATH

echo "Now running Airfoil OP2"
for i in $(seq 1 ${COUNT})
do
        echo "--------------------"
        airfoil/airfoil_plain/airfoil_dp_openmp new_grid.dat
        echo "--------------------"
done

echo "Now running Aero OP2"
for i in $(seq 1 ${COUNT})
do
        echo "--------------------"
        aero/aero_plain/aero_dp_openmp FE_grid.dat
        echo "--------------------"
done

cd ../../../../
