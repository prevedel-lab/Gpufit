To compile the project you can run:
```bash
cmake . -B build -G "Visual Studio 17 2022" -A x64 -DUSE_CUBLAS=True
cmake --build build --config Release --target ALL_BUILD --verbose
```