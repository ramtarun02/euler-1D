# prepare the compilation (ensure that we have a clean build/ directory)
Remove-Item .\build -Force -Recurse
New-Item -Name "build" -ItemType "directory"

# compile the solver into an object file
cl.exe /nologo /EHsc /std:c++17 /I. /c /O2 .\euler.cpp /Fo".\build\euler.obj"

# link object file into an executable
cl.exe /nologo .\build\euler.obj /Fe".\build\euler.exe"

# run the solver
Set-Location .\build
.\euler.exe
Set-Location ..