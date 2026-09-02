# Transport Calculation API

## How to Use
1. Enter the directory path containing the POSCAR files.
2. Enter the directory path containing the all the atomic pseudopotentials
3. Enter the directory path containing the hdf5-1.14.6 library.
4. Enter the file path of the INCAR file to use. Ensure the parameters are set to relaxation settings.
5. Enter the file path of the VASP binary file.
6. Enter the maximum calculation time (in minutes) for the relaxation, electron-phonon, and transport calculations.
7. Review the directories that the program will create and ensure they do not exist already.

## Keep Track of Calculation Progress
The program will create 3 log files: relaxed.log, phelel.log, and transport.log. Each log will write out the directories in which calculation will/have occurred. The letter in front of the directory is the status of the calculation.
```
  L: Loading in the calculation to be submitted, and is now in progress.
  P: The calculation passed.
  F: The calculation failed.
```
## NOTICE
**Use the absolute path for all file/directory paths for the best results.** \
The transport.log file is within the poscar_phelel folder and not in its own to reuse folders and reduce space.
