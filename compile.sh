g++ *.cpp -fdiagnostics-color=always -g ${file} -o robust_energy_path -std=c++17 -Wshadow -Wall -Wno-unused-result -DIL_STD -I/opt/ibm/ILOG/CPLEX_Studio_Community201/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio_Community201/concert/include -L/opt/ibm/ILOG/CPLEX_Studio_Community201/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio_Community201/concert/lib/x86-64_linux/static_pic -lilocplex -lconcert -lcplex -lm -lpthread -ldl -lemon