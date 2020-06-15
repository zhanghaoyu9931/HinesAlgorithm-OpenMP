if not exist ..\src\serial (
    echo "compile serial program..."
    g++ -std=c++11 -fopenmp ..\src\serial.cpp -o ..\src\serial
)
for /l %%i in (1,1,12) do ..\src\serial ..\data\case%%i.txt ..\sresult\res%%i.txt
pause