if not exist ..\src\parallel (
    echo "compile parallel program..."
    g++ -std=c++11 -fopenmp ..\src\parallel.cpp -o ..\src\parallel
)
for /l %%i in (1,1,12) do ..\src\parallel ..\data\case%%i.txt ..\presult\res%%i.txt
pause