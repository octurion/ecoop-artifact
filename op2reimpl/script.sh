echo "Now running Aero AoS"
echo "--------------------"
build/aero_aos FE_grid.dat /dev/null
echo "--------------------"

echo "Now running Aero partial SoA"
echo "--------------------"
build/aero FE_grid.dat /dev/null
echo "--------------------"

echo "Now running Aero partial SoA with opts"
echo "--------------------"
build/aero_v2 FE_grid.dat /dev/null
echo "--------------------"

echo "Now running Airfoil AoS"
echo "--------------------"
build/airfoil_aos new_grid.dat /dev/null
echo "--------------------"

echo "Now running Airfoil partial SoA"
echo "--------------------"
build/airfoil new_grid.dat /dev/null
echo "--------------------"

echo "Now running airfoil full SoA no opts"
echo "--------------------"
build/airfoil_v2 new_grid.dat /dev/null
echo "--------------------"

echo "Now running airfoil full SoA with opts"
echo "--------------------"
build/airfoil_v3 new_grid.dat /dev/null
echo "--------------------"
