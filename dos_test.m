for n = 1:10
    Str = sprintf('java -jar D:\\Download\\MATLAB\\phylonet_v2_4\\phylonet_v2_4\\phylonet_v2_4.jar rf -m D:\\Download\\MATLAB\\EvolDist\\Data\\Simulation_%d.dnd -e D:\\Download\\MATLAB\\EvolDist\\Data\\Simulation_ClustalX.ph -o D:\\Download\\MATLAB\\EvolDist\\Data\\RFdistance\\Simulation_%d.txt',n,n);
    [status, result] = dos(Str);
    disp(Str);
end