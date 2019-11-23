#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>

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
	std::cout << "Current path is " << std::filesystem::current_path() << std::endl << std::endl;

	using namespace SudokuSolving;
	//for (int i = 0; i < 16; i++)
	//{
	//	for (int j = 0; j < 16; j++)
	//	{
	//		std::cout << setColor((color)i, (color)j) << "hallo";
	//	}
	//	std::cout << std::endl;
	//}
	 
	while (true)
	{
		// Scan for sudoku files:
		std::vector<std::string> files;
		for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path()))
		{
			auto filename = entry.path().string();
			if (filename.rfind(".txt") != -1)
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
		if(s.fromFile(filename))
			continue;

		// Print it:
		s.out();
		initialCalcHints(&s);
		s.out_hints();

		// Load solvers:
		NakedSolver  nakedsol;
		HiddenSolver hiddensol;
		FishSolver	 fishsol;

		// Interactively solve the Sudoku:
		std::string input;
		while (!s.isSolved())
		{
			std::cout << "Choose solving technique... (n = naked, h = hidden, f = fish)" << std::endl;
			std::getline(std::cin, input);
			if (input == "")
				continue;
			int order = 3;
			if (input.size() > 1 && std::isdigit(input[1]))
				order = input[1] - '0';
			std::cout << "maximal order: " << order << std::endl;

			int v0 = s.getNumValues();
			int h0 = s.getNumHints();
			if (input[0] == 'n')
			{
				nakedsol.apply(&s, order);
			}
			else if (input[0] == 'h')
			{
				hiddensol.apply(&s, order);
			}
			else if (input[0] == 'f')
			{
				fishsol.apply(&s, order);
			}
			else
				continue;

			if(s.getNumValues() - v0 > 0)
				s.out();
			if(h0 - s.getNumHints() > 0 && !s.isSolved())
				s.out_hints();

			std::cout << "Filled " << s.getNumValues() - v0 << " cells." << std::endl;
			std::cout << "Eliminated " << h0 - s.getNumHints() << " hints." << std::endl;
		}
		std::cout << "Sudoku is solved!" << std::endl << std::endl;
		std::cin.get();
	}
	return 0;
}
