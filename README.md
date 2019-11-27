# SudokuSolver
An interactive step-by-step solver for regular Sudokus.

The executable reads regular (9x9) Sudokus from text files and offers a collection of algorithms to solve them. 
There is no GUI, only console with colored output. 

This little project is inspired by the great explanations on http://hodoku.sourceforge.net/en/techniques.php. 
The naming convention used in the console output was adapted from there and most of the included test sudokus are taken 
from the 'Solving techniques' section.   

### Supported solving techniques
* Locked Candidates
* Hidden Subsets, including Hidden Single
* Naked Subsets, including Last Digit
* Fish
  * Basic Fish: X-Wing (order 2), Swordfish (order 3), Jellyfish (order 4)
  * Exo Fins
  * Franken Fish and Mutant Fish 
* Brute Force

#### Planned solving techniques
* Locked Naked Subsets
* Chains
* Endo Fins
* Wings
* Uniqueness
* Siamese Fish

### Usage
Place Sudoku files in the same folder. After loading it, type one of the printed characters. You can append a number to specify the maximal order, the pattern is searched for. Example: h6 looks for all hidden subsets with at most 6 cells.

### Gallery
<p align="middle">
  <img width=450 src=/pics/naked_subsets.png>
  <img width=400 src=/pics/mutant_jellyfish.png>
</p>

### Compile
The code is pure C++, only using the standard library. Requires C++17 and std::map<>::contains from C++20.
Tested for Windows 10 and Ubuntu. The appearence of formatted console output may vary with the operating system.
