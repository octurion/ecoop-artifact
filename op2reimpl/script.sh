COUNT=20

echo "Now running Aero AoS"
for i in $(seq 1 ${COUNT})
do
    echo "--------------------"
    build/op2reimpl/aero_aos build/op2reimpl/FE_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Aero partial SoA"
for i in $(seq 1 ${COUNT})
do
    echo "--------------------"
    build/op2reimpl/aero build/op2reimpl/FE_grid.dat /dev/null
    echo "--------------------"
done

# echo "Now running Aero partial SoA with opts"
# echo "--------------------"
# build/aero_v2 FE_grid.dat /dev/null
# echo "--------------------"

echo "Now running Airfoil AoS"
for i in $(seq 1 ${COUNT})
do
    echo "--------------------"
    build/op2reimpl/airfoil_aos build/op2reimpl/new_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Airfoil partial SoA"
for i in $(seq 1 ${COUNT})
do
    echo "--------------------"
    build/op2reimpl/airfoil build/op2reimpl/new_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Airfoil full SoA"
for i in $(seq 1 ${COUNT})
do
    echo "--------------------"
    build/op2reimpl/airfoil_soa build/op2reimpl/new_grid.dat /dev/null
    echo "--------------------"
done
