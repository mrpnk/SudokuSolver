#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>

#define NOMINMAX 
#include "windows.h"
#ifdef _WINDOWS_
	
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
#endif

enum color
{
	black = 0,
	blue,
	green,
	lblue,
	red,
	purple,
	orange,
	lgray,
	gray,
	mblue,
	lgreen,
	turque,
	lred,
	lpurple,
	yellow,
	white
};
struct colorSetter { color fg; color bg; };
colorSetter setColor(color fg = white, color bg = black) { return colorSetter{ fg,bg }; }

std::ostream& operator<<(std::ostream& out, const colorSetter& f)
{
#ifdef _WINDOWS_
	SetConsoleTextAttribute(hConsole, f.fg + f.bg * 16);
#endif
	return out;
}

#include "sudoku.h"

bool is_number(const std::string& s)
{
	return !s.empty() && std::find_if(s.begin(),
		s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

int main(int argc, char* argv[])
{
	using namespace sudoku;

#if defined(_DEBUG)
	std::cout << "Current path is " << std::filesystem::current_path() << std::endl << std::endl;
#endif
 
	// Load solvers:
	LockedSolver     lockedsol; 
	HiddenSolver     hiddensol; 
	NakedSolver      nakedsol;
	FishSolver       fishsol;
	UniqueSolver     unisol;
	BruteforceSolver brutsol;

	while (true)
	{
		// Scan for Sudoku files:
		std::vector<std::string> files;
		for (auto& entry : std::filesystem::recursive_directory_iterator(std::filesystem::current_path()))
		{
			auto filename = entry.path().string();
			if (filename.rfind(".sudoku") != -1)
			{
				std::cout << files.size() << "   " << entry.path().filename() << std::endl;
				files.push_back(filename);
			}
		}

		// Select file:
		std::string filename;
		while (true)
		{
			std::cout << "Choose file..." << std::endl;
			std::string s;
			std::getline(std::cin, s);
			if (is_number(s))
			{
				int ix = std::stoi(s);
				if (ix >= 0 && ix < files.size())
				{
					filename = files[ix];
					break;
				}
				else std::cout << "No. This number is out of range!" << std::endl;
			}
			else
				std::cout << "No. Enter one of the numbers!" << std::endl;
		}

		// Load Sudoku:
		Sudoku s;
		if(s.loadFromFile(filename))
			continue;

		// Print it:
		std::cout<<s.valueRep<<std::endl;
		initialCalcCandidates(&s);
		std::cout<<s.candiRep<<std::endl;

		// Interactively solve the Sudoku:
		std::string input;
		while (!s.isSolved())
		{
			std::cout << "Choose solving technique..." << std::endl;
			std::getline(std::cin, input);
			if (input == "")
				continue;
			int order = 3;
			if (input.size() > 1 && std::isdigit(input[1]))
				order = input[1] - '0';
			std::cout << "maximal order: " << order << std::endl;

			int v0 = s.getNumValues();
			int h0 = s.getNumCandidates();
			if (input[0] == 'l')
			{
				lockedsol.apply(&s, order);
			}
			else if (input[0] == 'h')
			{
				hiddensol.apply(&s, order);
			}
			else if (input[0] == 'n')
			{
				nakedsol.apply(&s, order);
			}
			else if (input[0] == 'f')
			{
				fishsol.apply(&s, order);
			}
			else if (input[0] == 'u')
			{
				unisol.apply(&s, order);
			}
			else if (input[0] == 'b')
			{
				brutsol.apply(&s, order);
			}
			else if (input[0] == 'c')
			{
				break;
			}
			else if (input[0] == 'q')
			{
				return 0;
			}
			else
			{
				std::cout << "(l = locked, h = hidden, n = naked, f = fish, u = uniqueness, b = bruteforce; c = cancel, q = quit)" << std::endl;
				continue;
			}

			if(s.getNumValues() - v0 > 0)
				std::cout<<s.valueRep<<std::endl;
			if(h0 - s.getNumCandidates() > 0 && !s.isSolved())
				std::cout<<s.candiRep<<std::endl;

			std::cout << "Filled " << s.getNumValues() - v0 << " cells." << std::endl;
			std::cout << "Eliminated " << h0 - s.getNumCandidates() << " candidates." << std::endl;
		}
		if(s.isSolved())
			std::cout << "Sudoku is solved!" << std::endl << std::endl;
	}
	std::cout<<"Goodbye..."<<std::endl;
	return 0;
}
