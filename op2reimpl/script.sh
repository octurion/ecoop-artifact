echo "CPU information:"
echo "--------------------"
cat /proc/cpuinfo
echo "--------------------"

echo "Now running Aero AoS"
for i in {0..19}
do
    echo "--------------------"
    build/aero_aos FE_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Aero partial SoA"
for i in {0..19}
do
    echo "--------------------"
    build/aero FE_grid.dat /dev/null
    echo "--------------------"
done

# echo "Now running Aero partial SoA with opts"
# echo "--------------------"
# build/aero_v2 FE_grid.dat /dev/null
# echo "--------------------"

echo "Now running Airfoil AoS"
for i in {0..19}
do
    echo "--------------------"
    build/airfoil_aos new_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Airfoil partial SoA"
for i in {0..19}
do
    echo "--------------------"
    build/airfoil new_grid.dat /dev/null
    echo "--------------------"
done

echo "Now running Airfoil full SoA"
for i in {0..19}
do
    echo "--------------------"
    build/airfoil_soa new_grid.dat /dev/null
    echo "--------------------"
done

# echo "Now running airfoil full SoA no opts"
# echo "--------------------"
# build/airfoil_v2 new_grid.dat /dev/null
# echo "--------------------"

# echo "Now running airfoil full SoA with opts"
# echo "--------------------"
# build/airfoil_v3 new_grid.dat /dev/null
# echo "--------------------"
