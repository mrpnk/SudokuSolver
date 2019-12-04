#pragma once
#include <vector>
#include <deque>
#include <tuple>
#include <map>
#include <set>
#include <array>
#include <assert.h>

namespace SudokuSolving
{
	short count_bit_o(unsigned int i)
	{
		// Java: use >>> instead of >>
		// C or C++: use uint32_t
		i = i - ((i >> 1) & 0x55555555);
		i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
		return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}
	unsigned short bitcount(unsigned short x)
	{
		x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
		x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
		x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
		//x = (x & 0x0000FFFF) + ((x >> 16) & 0x0000FFFF);
		return x;
	}
	int get_last_set(uint16_t bits)
	{
		if (bits == 0)
			return -1;
		int value = 0;
		while ((bits & 1) == 0)
		{
			bits >>= 1;
			value++;
		}
		return value;
	}

	struct Cell
	{
		unsigned short value;
		unsigned short hints;
	};
	class Sudoku
	{
		int n_;
		int nn_;
	public:
		struct Indexer
		{
			friend class Sudoku;
			enum Type
			{
				row,
				clm,
				box
			};
			virtual int operator ()(int i) = 0;
			virtual Type type() = 0;
			int j()
			{
				return j_;
			}
			Indexer(int j, int n, int nn) : j_{ j }, n{ n }, nn{ nn }{}
		protected:
			int j_;
			int n;
			int nn;
		};

		std::vector<Cell> cells; // rowwise, as always...
		Indexer::Type types[3] = { Indexer::Type::row, Indexer::Type::clm, Indexer::Type::box };

		struct Marker
		{
			int markGroup;
		};
		std::map<std::tuple<short, short, short>, Marker> markers; // x,y,hint(0=value)

		Sudoku(int n = 0)
		{
			init(n);
		}
		void init(int n)
		{
			n_ = n;
			nn_ = n * n;

			for (int j = 0; j < nn_; j++)
			{
				row_indexer.push_back({ j, n_, nn_ });
				clm_indexer.push_back({ j, n_, nn_ });
				box_indexer.push_back({ j, n_, nn_ });
			}

			cells.resize(nn_ * nn_);
			for (int i = 0; i < nn_ * nn_; i++)
			{
				cells[i].value = 0;
				cells[i].hints = 0xffff >> (16 - nn_);
			}
		}

		int n()
		{
			return n_;
		}
		int nn()
		{
			return nn_;
		}
		int nnn()
		{
			return nn_*n_;
		}
		void set(int x, int y, int value)
		{
			cells[y * nn_ + x].value = 1 << (value - 1);
			cells[y * nn_ + x].hints = 0;
		}
		int get(int x, int y)
		{
			unsigned short bits = cells[y * nn_ + x].value;
			return get_last_set(bits) + 1;
		}
		bool hint(int x, int y, int hint)
		{
			return (bool)(cells[y * nn_ + x].hints & (1 << (hint - 1)));
		}
		void out()
		{
			using namespace std;
			std::vector<color> cols = { turque,orange,lblue,lpurple,green,red,yellow,purple,blue,gray,gray };
			int maxiter = 0;
			for (int y = 0; y < nn_ + 1; y++)
			{
				cout << "   ";
				if (y % n_ == 0)
					cout << string(nn_ * 2 + n_ + (n_ + 1), '-') << endl << "   ";
				if (y < nn_)
				{
					for (int x = 0; x < nn_; x++)
					{
						if (x % n_ == 0)
							cout << "| ";
						int value = get(x, y);
						bool mark = markers.contains({ x,y,0 });
						int iter = markers[{x, y, 0}].markGroup;
						maxiter = max(maxiter, iter);
						if (value != 0)
							cout << setColor(mark ? cols[iter] : white) << value << setColor(white) << " ";
						else
							cout << "  ";
					}
					cout << "|" << endl;
				}
			}
			cout << endl << "Iterations: ";
			for (int i = 0; i <= maxiter; i++)
			{
				cout << setColor(cols[i]) << i << setColor() << (i == maxiter ? "." : ", ");
			}
			cout << endl << endl;
		}
		void out_hints()
		{
			std::vector<color> foreColors = { turque, orange, red, lgreen };
			std::vector<color> backColors = { black, gray, blue, purple };
			int rowhl = 0, colhl = 0, boxhl = 0;

			using namespace std;
			bool marked = false;
			for (int y = 0; y < nn_ + 1; y++)
			{
				if (y % n_ == 0)
					cout << string(nn_ * (n_ + 2) + n_ * 2 + (n_ + 1), '-') << endl;
				rowhl = markers.contains({ -1,y,0 }) ? 1 : 0;
				if (y < nn_)
				{
					for (int iy = 0; iy < n_; iy++)
					{
						for (int x = 0; x < nn_; x++)
						{
							colhl = markers.contains({ x,-1,0 }) ? 1 : 0;
							boxhl = markers.contains({ -1,-1,getBox(nn_*y+x) }) ? 1 : 0;
							if (x == 0)
								std::cout << setColor() << "|" << setColor(white, backColors[rowhl]) << "  ";
							else if (x % n_ == 0)
								std::cout << setColor(white, backColors[rowhl]) << "|  ";

							for (int ix = 0; ix < n_; ix++)
							{
								short h = iy * n_ + ix + 1;
								bool mark = markers.contains({ x,y,h });
								marked |= mark;
								auto& m = markers[{x, y, h}];
								bool forced = !hint(x, y, h) && mark;
								std::cout << setColor(mark ? foreColors[m.markGroup] : white, forced ? yellow : backColors[rowhl + colhl + boxhl]) <<
									((hint(x, y, h) || mark) ? to_string(h) : " ") << setColor(white, backColors[rowhl]);
							}

							std::cout << "  ";
						}
						cout << setColor() << "|" << endl;
					}

					if ((y + 1) % n_ != 0)
						cout << endl;
				}
			}
			if (marked)
				cout << "(" << setColor(foreColors[0]) << "base" << setColor() << ", " << setColor(foreColors[1]) << "cover" << setColor() << ", " << setColor(foreColors[2]) << "eliminated" << setColor() << ", " << setColor(foreColors[3]) << "fin" << setColor() << ")" << endl;
			std::cout << endl;
		}

		Indexer* getIndexer(Indexer::Type t, int j)
		{
			if (t == Indexer::Type::row)
				return &row_indexer[j];
			else if (t == Indexer::Type::clm)
				return &clm_indexer[j];
			else if (t == Indexer::Type::box)
				return &box_indexer[j];
			else return nullptr;
		}

		bool fromFile(std::string filename)
		{
			try
			{
				std::ifstream infile(filename);

				std::string line;
				std::getline(infile, line);
				init(std::stoi(line));

				int y = 0;
				while (std::getline(infile, line))
				{
					int index = 0;
					int x = 0;
					while (index < line.size()) // TODO multiple character values, spaces
					{
						if (line[index] == ',')
							x++;
						else
						{
							int value = std::stoi(line.substr(index, 1));
							set(x, y, value);
						}
						index++;
					}
					y++;
					if (y >= nn())
						break;
				}
				return false;
			}
			catch (std::exception & e)
			{
				std::cout << setColor(red) << "Parser error in sudoku file " << filename << ". Discription: " << e.what() << setColor() << std::endl;
				return true;
			}
		}

		bool isSolved()
		{
			for (auto& a : cells)
				if (a.value == 0)
					return false;
			return true;
		}

		int getNumValues()
		{
			int ret = 0;
			for (auto& a : cells)
				ret += (a.value != 0);
			return ret;
		}
		int getNumHints()
		{
			int ret = 0;
			for (auto& a : cells)
				ret += bitcount(a.hints);
			return ret;
		}

		int getRow(int cell)
		{
			return cell / (nn_); // get row index from cell index
		}
		int getClm(int cell)
		{
			return cell % (nn_); // get column index from cell index
		}
		int getBox(int cell)
		{
			return (cell / n_) % n_ + n_ * (cell / (nn_*n_)); // get box index from cell index
		}
		int getHouse(Indexer::Type t, int cell) 
		{
			switch (t)
			{
			case Indexer::row:
				return getRow(cell);
			case Indexer::clm:
				return getClm(cell);
			case Indexer::box:
				return getBox(cell);
			}
		}
	private:
		struct rowIndexer : public Indexer {
			using Indexer::Indexer;
			int operator()(int i) override
			{
				return i + j_ * nn;
			}
			Type type() override { return Indexer::Type::row; }
		};
		std::vector<rowIndexer> row_indexer;

		struct clmIndexer : public Indexer {
			using Indexer::Indexer;
			int operator()(int i) override
			{
				return i * nn + j_;
			}
			Type type() override { return Indexer::Type::clm; }
		};
		std::vector<clmIndexer> clm_indexer;

		struct boxIndexer: public Indexer {
			using Indexer::Indexer;
			int operator()(int i) override
			{
				int base = (j_ % n) * n + (j_ / n) * n * nn; // base aus box-index j
				int box = base + i % n + (i / n) * nn;
				return box;
			}
			Type type() override { return Indexer::Type::box; }
		};
		std::vector<boxIndexer> box_indexer;
	};

	std::string descr(Sudoku::Indexer::Type t)
	{
		return t == Sudoku::Indexer::Type::row ? "r" : t == Sudoku::Indexer::Type::clm ? "c" : "b";
	}
	std::string getTextFromBits(uint16_t bits, int baseIndex, std::string seperator = ", ")
	{
		std::string ret;
		for (int value = 0; value < 16; value++)
		{
			if ((bits & (1 << value)) != 0)
				ret += (ret == "" ? "" : seperator) + std::to_string(value + baseIndex);
		}
		return ret;
	}

	int eliminateHints(std::vector<Cell>& cells, int nn, uint16_t hints, Sudoku::Indexer* indexer0,
		Sudoku::Indexer* indexer1 = nullptr, Sudoku::Indexer* indexer2 = nullptr,
		uint16_t positions = 0xffff)
	{
		int numChanges = 0;
		for (int i = 0; i < nn; i++)
		{
			if (((1 << i) & positions) == 0)
				continue;

			// eliminate hints in first indexer (row propably)
			if (indexer0)
			{
				uint16_t& cellhints = cells[(*indexer0)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}

			// eliminate hints in second indexer (column propably)
			if (indexer1)
			{
				uint16_t& cellhints = cells[(*indexer1)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}

			// eliminate hints in third indexer (box propably)
			if (indexer2)
			{
				uint16_t& cellhints = cells[(*indexer2)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}
		}
		return numChanges;
	}
	int eliminateHints(Sudoku* s, uint16_t hints, Sudoku::Indexer* indexer0,
		Sudoku::Indexer* indexer1 = nullptr, Sudoku::Indexer* indexer2 = nullptr,
		uint16_t positions = 0xffff)
	{
		int numChanges = 0;
		for (int i = 0; i < s->nn(); i++)
		{
			if (((1 << i) & positions) == 0)
				continue;

			// eliminate hints in first indexer (row propably)
			if (indexer0)
			{
				uint16_t& cellhints = s->cells[(*indexer0)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}

			// eliminate hints in second indexer (column propably)
			if (indexer1)
			{
				uint16_t& cellhints = s->cells[(*indexer1)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}

			// eliminate hints in third indexer (box propably)
			if (indexer2)
			{
				uint16_t& cellhints = s->cells[(*indexer2)(i)].hints;
				if ((cellhints & hints) != 0)
					numChanges++;
				cellhints &= ~hints;
			}
		}
		return numChanges;
	}
	void initialCalcHints(struct Sudoku* s)
	{
		for (int cell = 0; cell < s->nn() * s->nn(); cell++)
		{
			if (s->cells[cell].value != 0)
			{
				eliminateHints(s, s->cells[cell].value,
					s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
					s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
					s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell)));
			}
		}
	}

	bool intersects(Sudoku* s, Sudoku::Indexer* indexer0, Sudoku::Indexer* indexer1, uint16_t hints)
	{
		// check if any of the given hints are in intersections between structures
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			int index0 = (*indexer0)(i0);
			for (int i1 = 0; i1 < s->nn(); i1++)
			{
				if (index0 == (*indexer1)(i1))
				{
					if ((s->cells[index0].hints & hints) != 0x0)
						return true;
				}
			}
		}
		return false;
	}
	bool intersects(Sudoku* s, Sudoku::Indexer* indexer0, std::vector<bool> positions, uint16_t hints)
	{
		// check if any of the given hints are in intersections between structures
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			int index0 = (*indexer0)(i0);
			if (positions[index0])
			{
				if ((s->cells[index0].hints & hints) != 0x0)
					return true;
			}

		}
		return false;
	}
	bool intersects(std::vector<bool> positions0, std::vector<bool> positions1)
	{
		// check if any of the given hints are in intersection
		for (int cell = 0; cell < positions0.size(); cell++)
		{
			if (positions0[cell] && positions1[cell])
			{
				return true;
			}

		}
		return false;
	}

	std::vector<bool> merge(Sudoku* s, Sudoku::Indexer* indexer0, std::vector<bool> positions, uint16_t hints)
	{
		std::vector<bool> ret = positions;

		// check if any of the given hints are in intersections between structures
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			int index0 = (*indexer0)(i0);
			if ((s->cells[index0].hints & hints) != 0x0)
			{
				ret[index0] = true;
			}
		}
		return ret;
	}

	bool isempty(std::vector<bool>& positions)
	{
		for (bool b : positions)
			if (b)
				return false;
		return true;
	}
	int count(const std::vector<bool>& positions)
	{
		int count{ 0 };
		for (const auto& a : positions)
			if (a) count++;

		return count;
	}

	void visible(Sudoku* s, const std::vector<bool>& cells, std::vector<bool>& positions)
	{
		// unsets all positions that are not visible from any of the given cells
		for (int j = 0; j < s->nn() * s->nn(); j++)
		{
			if (!cells[j])
				continue;
			int cr = s->getRow(j);
			int cc = s->getClm(j);
			int cb = s->getBox(j);
			for (int i = 0; i < s->nn() * s->nn(); i++)
			{
				auto temp = positions[i];
				positions[i] = positions[i] && ((cr == s->getRow(i)) || (cc == s->getClm(i)) || (cb == s->getBox(i)));
				if (temp != positions[i])
					int idsf = 678;
			}
		}
	}
	void visible(Sudoku* s, const std::set<short>& cells, std::set<short>& positions, bool includeCells)
	{
		assert(cells.size() <= 16);
		// sets positions to the cells that are visible by all given cells
		std::map<short, short> os;
		int j = 0;
		for (auto c : cells)
		{
			for (auto t : s->types)
			{
				auto indexer = s->getIndexer(t, s->getHouse(t, c));
				for (int i = 0; i < s->nn(); i++)
					os[(*indexer)(i)] |= (1 << j);
			}
			j++;
		}
		for (auto& a : os)
			if (bitcount(a.second) == cells.size() && (!cells.contains(a.first) || includeCells))
				positions.insert(a.first);
	}

	bool subtract(const std::vector<bool>& positions0, const std::vector<bool>& positions1, std::vector<bool>& only0, std::vector<bool>& only1, bool breakOnOnly1) // returns if 1 is entirely in 0
	{
		for (int cell = 0; cell < positions0.size(); cell++)
		{
			if (positions0[cell] && !positions1[cell])
			{
				only0[cell] = true;
			}
			else if (!positions0[cell] && positions1[cell])
			{
				only1[cell] = true;
				if (breakOnOnly1)
					return false;
			}
		}
		return true;
	}
	void subtract(Sudoku* s, Sudoku::Indexer* indexer0, Sudoku::Indexer* indexer1, uint16_t& positions0, uint16_t hintmask = 0xffff)
	{
		// sets excactly those positions0 that are not part of indexer1 but have at least one given hint.
		positions0 = 0x0;
		for (int i = 0; i < s->nn(); i++)
		{
			auto cell = (*indexer0)(i);
			if (s->cells[cell].hints & hintmask)
			{
				uint16_t p = (1u << i);
				positions0 |= p;
				for (int j = 0; j < s->nn(); j++)
				{
					if ((*indexer1)(j) == cell)
					{
						positions0 &= ~p;
						break;
					}
				}
			}
		}
	}
	uint16_t containingAll(Sudoku* s, Sudoku::Indexer* indexer, uint16_t hints)
	{
		// get the subindices of those cells, that contain all given hints
		uint16_t positions = 0x0;
		for (int i = 0; i < s->nn(); i++)
		{
			uint16_t temp = s->cells[(*indexer)(i)].hints;
			if ((temp & hints) == hints)
				positions |= (1 << i);
		}
		return positions;
	}
	uint16_t containingAny(Sudoku* s, Sudoku::Indexer* indexer, uint16_t hints)
	{
		// get the subindices of those cells, that contain any of the given hints
		uint16_t positions = 0x0;
		for (int i = 0; i < s->nn(); i++)
		{
			uint16_t temp = s->cells[(*indexer)(i)].hints;
			if (temp & hints)
				positions |= (1 << i);
		}
		return positions;
	}
	Sudoku::Indexer* secondHouse(Sudoku* s, Sudoku::Indexer* indexer, uint16_t positions)
	{
		// returns the indexer of the second house, all the given positions of the given indexer are in, if there is one.
		uint16_t secondHouse1 = 0xffff, secondHouse2 = 0xffff;
		for (int i = 0; i < s->nn(); i++)
		{
			if ((positions & (1 << i)) != 0)
			{
				switch (indexer->type())
				{
					case Sudoku::Indexer::row:
					case Sudoku::Indexer::clm:
						secondHouse1 &= (1 << s->getBox((*indexer)(i)));
						secondHouse2 = 0x0;
						break;
					case Sudoku::Indexer::box:
						secondHouse1 &= (1 << s->getRow((*indexer)(i)));
						secondHouse2 &= (1 << s->getClm((*indexer)(i)));
						break;
				}
				
				if (!secondHouse1 && !secondHouse2)
					return nullptr;
			}
		}
		switch (indexer->type())
		{
		case Sudoku::Indexer::row:
		case Sudoku::Indexer::clm:
			return s->getIndexer(Sudoku::Indexer::box, get_last_set(secondHouse1));
		case Sudoku::Indexer::box:
			return secondHouse1 ? s->getIndexer(Sudoku::Indexer::row, get_last_set(secondHouse1)) : s->getIndexer(Sudoku::Indexer::clm, get_last_set(secondHouse2));;
		}
	}
	template<typename T>
	std::set<T> bitindex(T* t, unsigned short idx, int max = 16)
	{
		std::set<T> ret;
		for (int i = 0; i < max; i++)
		{
			if (idx & (1 << i))
				ret.insert(t[i]);
		}
		return ret;
	}

	unsigned short bitintersection(std::set<unsigned short> hs)
	{
		// returns the intersection of all given bitmasks
		unsigned short ret = 0xffff;
		for (unsigned short s : hs)
			ret &= s;
		return ret;
	}
	unsigned short bitunion(std::set<unsigned short> hs)
	{
		// returns the intersection of all given bitmasks
		unsigned short ret = 0x0;
		for (unsigned short s : hs)
			ret |= s;
		return ret;
	}

	class SudokuSolver
	{
	public:
		virtual void apply(Sudoku* s, int maxorder) const = 0;
	};

	class LockedSolver : public SudokuSolver
	{
		struct lockedResult
		{
			Sudoku::Indexer* indexer;
			uint16_t hints;
			uint16_t positions;
		};
	protected:
		void find_locked(Sudoku* s, Sudoku::Indexer* indexer, std::vector<lockedResult>& results) const
		{
			for (short h = 0; h < s->nn(); h++)
			{
				uint16_t thisPositions = containingAll(s, indexer, (1 << h));
				if (thisPositions == 0x0)
					continue;

				// watch out for a second house all occurences of hint h are contained in.
				if (auto indexer2 = secondHouse(s, indexer, thisPositions))
				{
					uint16_t locked_positions;
					subtract(s, indexer2, indexer, locked_positions, 1u << h);
					if(locked_positions)
						results.push_back({ indexer2, 1u << h, locked_positions });
				}
			}
		}
		void handle_locked(Sudoku* s, const lockedResult& r, bool& changed, int iter) const
		{
			// erase every other hints from these positions
			if ((r.positions & (r.positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*r.indexer)(get_last_set(r.positions));
				s->cells[cell].value = r.hints;
				s->cells[cell].hints = 0x0;
				eliminateHints(s, r.hints,
					s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
					s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
					s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell)));
				changed = true;
				std::cout << " -> set " << getTextFromBits(r.hints, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateHints(s, r.hints, r.indexer, nullptr, nullptr, r.positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there z hints in a structure that distribute over m <= z cells?
			s->markers.clear();
			std::vector<lockedResult> results;
			bool changed = true;
			int iter = 1;
			while (changed)
			{
				changed = false;
				for (Sudoku::Indexer::Type t : s->types)
				{
					for (int j = 0; j < s->nn(); j++)
					{
						// look for hidden z-tuples in j'th structure of type t
						Sudoku::Indexer* indexer = s->getIndexer(t, j);
						find_locked(s, indexer, results);
						for (auto& r : results)
						{
							std::cout << "found Locked Candidates in " << descr(t) << j + 1 << "/" << descr(r.indexer->type()) << r.indexer->j() + 1
								<< ": hints " << getTextFromBits(r.hints, 1) << " on pos " << getTextFromBits(r.positions, 1);
							handle_locked(s, r, changed, iter);
							std::cout << std::endl;
						}
						results.clear();
					}
				}
				iter++;
			}
		}

	};

	class HiddenSolver : public SudokuSolver
	{
	protected:
		void find_hidden(Sudoku* s, Sudoku::Indexer* indexer, int z, int m, std::vector<uint16_t>& results,
			uint16_t hints = 0x0, uint16_t positions = 0x0, int start_hint = 0) const
		{
			if (z == 0)
			{
				// check if there are candidates to delete:
				uint16_t delcands = 0x0;
				for (int i = 0; i < s->nn(); i++)
					if (positions & (1 << i))
						delcands |= s->cells[(*indexer)(i)].hints;
				delcands &= ~hints;
				if (m == 1 || (bool)delcands)
				{
					results.push_back(hints);
					results.push_back(positions);
				}
				return;
			}

			for (short h = start_hint; h < s->nn(); h++)
			{
				uint16_t thisPositions = containingAll(s, indexer, (1 << h));
				if (thisPositions == 0x0)
					continue;

				// merge positions in house of given indexer
				uint16_t newpositions = positions | thisPositions;
				auto nposis = bitcount(newpositions);
		
				if (nposis > m)
					continue;

				find_hidden(s, indexer, z - 1, m, results, hints | (1 << h), newpositions, h + 1);
			}
		}
		void handle_hidden(Sudoku* s, Sudoku::Indexer* indexer, uint16_t positions, uint16_t hints, bool& changed, int iter) const
		{
			// erase every other hints from these positions
			if ((positions & (positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*indexer)(get_last_set(positions));
				s->cells[cell].value = hints;
				s->cells[cell].hints = 0x0;
				eliminateHints(s, hints,
					s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
					s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
					s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell)));
				changed = true;
				std::cout << " -> set " << getTextFromBits(hints, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateHints(s, ~hints, indexer, nullptr, nullptr, positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there z hints in a structure that distribute over m <= z cells?
			std::string subsetNames[] = { "Single", "Pair", "Triple", "Quadrupel", "Quintuple", "Sextuple", "Septuple", "Octuple" };
			s->markers.clear();
			std::vector<uint16_t> results; // first the hints, second sub-index positions
			bool changed = true;
			int iter = 1;
			while (changed)
			{
				changed = false;
				for (int z = 1; z <= maxorder; z++)
				{
					for (Sudoku::Indexer::Type t : s->types)
					{
						for (int j = 0; j < s->nn(); j++)
						{
							// look for hidden z-tuples in j'th structure of type t
							Sudoku::Indexer* indexer = s->getIndexer(t, j);
							find_hidden(s, indexer, z, z, results);
							for (int r = 0; r < results.size() / 2; r += 2)
							{
								std::cout << "found Hidden " << subsetNames[z-1] << " in " << descr(t)
									<< j + 1 << ": hints " << getTextFromBits(results[r], 1) << " on pos " << getTextFromBits(results[r + 1], 1);
								handle_hidden(s, indexer, results[r + 1], results[r], changed, iter);
								std::cout << std::endl;
							}
							results.clear();
						}
					}
				}
				iter++;
			}
		}

	};

	class NakedSolver : public SudokuSolver
	{
	protected:
		void find_naked(Sudoku* s, Sudoku::Indexer* indexer,
			int z, int m, std::vector<uint16_t>& results,
			uint16_t hints = 0x0, uint16_t positions = 0x0, int start_j = 0) const
		{
			if (z == 0)
			{
				// check if there are candidates to delete:
				uint16_t delposis = containingAny(s, indexer, hints) & ~positions;
				if (m == 1 || (bool)delposis)
				{
					results.push_back(hints);
					results.push_back(positions);
				}
				return;
			}

			for (short j = start_j; j < s->nn(); j++)
			{
				if (s->cells[(*indexer)(j)].hints == 0x0)
					continue;

				// merge hints in structure of given indexer
				uint16_t newhints = hints | s->cells[(*indexer)(j)].hints;
				auto nhints = bitcount(newhints);
				if (nhints > m)
					continue;

				find_naked(s, indexer, z - 1, m, results, newhints, positions | (1 << j), j + 1);
			}
		}
		void handle_naked(Sudoku* s, Sudoku::Indexer* indexer,
			uint16_t positions, uint16_t hints, bool& changed, int iter) const
		{
			// erase those hints from elsewhere in the structure
			if ((positions & (positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*indexer)(get_last_set(positions));
				s->cells[cell].value = hints;
				eliminateHints(s, hints,
					s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
					s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
					s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell)));
				changed = true;
				std::cout << " -> set " << getTextFromBits(hints, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateHints(s, hints, indexer, nullptr, nullptr, ~positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there m cells in a structure that share z <= m hints? Then remove these hints from elsewhere in the structure (here: m = z)
			std::string subsetNames[] = { "Single", "Pair", "Triple", "Quadrupel", "Quintuple", "Sextuple", "Septuple", "Octuple" };
			s->markers.clear();
			std::vector<uint16_t> results; // first the hints, second sub-index positions
			bool changed = true;
			int iter = 1;
			while (changed)
			{
				changed = false;
				for (int z = 1; z <= maxorder; z++)
				{
					for (Sudoku::Indexer::Type t : s->types)
					{
						for (int j = 0; j < s->nn(); j++)
						{
							// look for naked z-tuples in j'th structure of type t
							Sudoku::Indexer* indexer = s->getIndexer(t, j);
							find_naked(s, indexer, z, z, results);
							for (int r = 0; r < results.size() / 2; r += 2)
							{
								std::cout << "found Naked " << subsetNames[z - 1] << " in " << descr(t)
									<< j + 1 << ": hints " << getTextFromBits(results[r], 1) << " on pos " << getTextFromBits(results[r + 1], 1);
								handle_naked(s, indexer, results[r + 1], results[r], changed, iter);
								std::cout << std::endl;
							}
							results.clear();
						}
					}
				}
				iter++;
			}
		}

	};

	class FishSolver : public SudokuSolver
	{
		struct fishResult
		{
			std::vector<bool> positions;
			int finns;
		};

	protected:
		mutable std::vector<bool> basePositions;

		void find_fish(Sudoku* s, uint16_t hints, int z, int m, std::vector<bool> positions,
			std::deque<std::pair<Sudoku::Indexer::Type, int>> sets, std::vector<fishResult>& results,
			std::vector<std::pair<Sudoku::Indexer::Type, int>>& locations, bool finned, bool franken, int& ncombis,
			const int maxnum = -1, bool cover = false, std::array<int, 3> start_j = { 0,0,0 }, uint32_t baseSets = 0x0) const
		{
			if (z == 0)
			{
				if (cover)
				{
					ncombis++;
					// do the cover sets cover all given hints in the base sets?
					auto onlyCover = std::vector<bool>(s->nn() * s->nn(), false);
					auto onlyBase = std::vector<bool>(s->nn() * s->nn(), false);
					if (subtract(positions, basePositions, onlyCover, onlyBase, !finned))
					{
						// Select only the cover hints, that are visible by all fins
						auto temp = onlyCover;
						visible(s, onlyBase, onlyCover);

						if (!isempty(onlyCover))
						{
							visible(s, onlyBase, temp);

							// mark the hints, override the old markers
							s->markers.clear();
							for (int i = 0; i < positions.size(); i++)
								if (positions[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(hints) + 1}] = { 1 }; // cover candidates
							for (int i = 0; i < basePositions.size(); i++)
								if (basePositions[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(hints) + 1}] = { 0 }; // base candidates
							for (int i = 0; i < onlyBase.size(); i++)
								if (onlyBase[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(hints) + 1}] = { 3 }; // fins
							for (int i = 0; i < onlyCover.size(); i++)
								if (onlyCover[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(hints) + 1}] = { 2 }; // deletions
							for (auto& se : sets)
								if (se.first == Sudoku::Indexer::row)
									s->markers[{-1, se.second, 0}] = { 0 };
								else if(se.first == Sudoku::Indexer::clm)
									s->markers[{se.second, -1, 0}] = { 0 };
								else
									s->markers[{-1, -1, se.second}] = { 0 };

							results.push_back({ onlyCover,count(onlyBase) });
							for (int zz = 0; zz < 2 * m; zz++)
							{
								locations.push_back(sets[zz]);
							}
						}
					}
					return;
				}
				else
				{
					// Reset the values. The base sets are found, Now look for cover sets.
					basePositions = positions;
					//basePositions.flip();
					positions = std::vector<bool>(s->nn() * s->nn(), false);

					z = m;
					cover = true;
					start_j[0] = start_j[1] = start_j[2] = 0;
				}
			}

			// for every structure
			for (Sudoku::Indexer::Type t : s->types)
			{
				if (!franken && t == Sudoku::Indexer::Type::box) // without franken, dont consider for boxes
					continue;
				//if (!cover && !sets.empty() && t != sets[0].first) // dont allow base sets to be of more than one type
				//	continue;
				//if (cover && t == sets[0].first) // dont allow base- and cover-sets to be of same type
				//	continue;
				for (int j = start_j[t]; j < s->nn(); j++)
				{
					if (cover && (baseSets & (1 << (j + t * s->nn()))))
						continue;

					// get indexer
					Sudoku::Indexer* indexer = s->getIndexer(t, j);

					// if there is no given hint in this structure or a hint is also in an already chosen structure, this is invalid
					uint16_t thisPositions = containingAll(s, indexer, hints);
					if (thisPositions == 0x0 || intersects(s, indexer, positions, hints))
						continue;

					// merge positions in structure of given indexer
					auto newpositions = merge(s, indexer, positions, hints);
					/*if (count(newpositions) > m)
						continue;*/

					auto newSets = sets;
					newSets.push_back({ t, j });
					start_j[t] = j;
					// go on to find next structure
					find_fish(s, hints, z - 1, m, newpositions, newSets, results, locations, finned, franken, ncombis, maxnum, cover, start_j, cover ? baseSets : baseSets | (1<<(j+t*s->nn())));
					newSets.pop_back();

					if (results.size() == maxnum)
						return;
				}
			}
		}

		void handle_fish(Sudoku* s, std::vector<bool> positions, uint16_t hints, bool& changed) const
		{
			// erase every other hints from these positions
			std::cout << " -> eliminate " << count(positions) << " occurences";
			for (int cell = 0; cell < s->nn() * s->nn(); cell++)
			{
				if (positions[cell])
				{
					// erase
					uint16_t& cellhints = s->cells[cell].hints;
					if ((cellhints & ~hints) != 0)
						changed = true;
					cellhints &= ~hints;

					// set
					if ((cellhints & (cellhints - 1)) == 0) /// only one bit set?
					{
						std::cout << " -> set " << getTextFromBits(cellhints, 1) << " at r" << s->getRow(cell) + 1 << "c" << s->getClm(cell) + 1;
						s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { 0 };

						s->cells[cell].value = cellhints;
						eliminateHints(s, cellhints,
							s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
							s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
							s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell)));
						changed = true;
					}
				}
			}
		}


	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			std::string fishnames[] = { "X-Wing", "Swordfish", "Jellyfish" };
			s->markers.clear();
			std::vector<fishResult> results;
			std::vector<std::pair<Sudoku::Indexer::Type, int>> locs;

			bool changed = true;
			while (changed)
			{
				changed = false;
				for (int z = 2; z <= maxorder; z++)
				{
					if (z > 4)
					{
						std::cout << "Maximal fish order is 4. For 9x9 Sudokus this covers all cases. See http://hodoku.sourceforge.net/de/tech_fishb.php#bf5" << std::endl;
						break;
					}
					int ncombis = 0;
					int searchMS = 0;
					for (short h = 0; h < s->nn(); h++)
					{
						auto startT = std::chrono::high_resolution_clock::now();
						find_fish(s, (1 << h), z, z, std::vector<bool>(s->nn() * s->nn(), false), {}, results, locs, true, true, ncombis, 1);
						searchMS += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startT).count();
						for (int r = 0; r < results.size(); r += 1)
						{
							std::cout << "Order " << z << " - search (until hit) took " << searchMS << " ms (" << ncombis << " combinations)" << std::endl;
							std::string attribute = (results[r].finns == 0 ? "" : std::to_string(results[r].finns) + "-finned "), sets;
							std::string oldtype, type;
							int category = 0; unsigned flags=0;
							for (int zz = 0; zz <= 2 * z; zz++)
							{
								if (zz % z == 0)
								{
									category = max(max(category, flags&4), 2*(flags % 4 == 3)); // try to distinguish franken from mutant fish
								}
								if (zz == z)
								{
									sets += "/";
									oldtype = "";
									flags = 0;
								}
								if (zz == 2 * z)
									break;
								auto t = locs[r * 2 * z + zz].first;
								flags |= (1u << t);
								type = descr(t);
								sets += (type != oldtype ? type : "") + std::to_string(locs[r * 2 * z + zz].second + 1);
								oldtype = type;
							}
							attribute = attribute + (category == 0 ? "" : category == 1 ? "Franken " : "Mutant ");
							std::cout << "found " << attribute << fishnames[z - 2] << " of hint " << (h + 1) << " in " << sets;
							handle_fish(s, results[r].positions, (1 << h), changed);
							std::cout << std::endl;
							if (changed)
								return;
						}
						results.clear();
						locs.clear();
					}
					std::cout << "Order " << z << " - search took " << searchMS << " ms  (" << ncombis << " combinations)" << std::endl;
				}
				break;
			}
		}
	};

	class UniqueSolver : public SudokuSolver
	{
		struct urResult
		{
			unsigned short urrows;
			unsigned short urclms;
			unsigned short urhints;
			int urtype;
			std::set<short> ercells;
			unsigned short erhints;
		};

		unsigned short additionals(unsigned short h[4], unsigned short& leastCommon) const
		{

			// finds all cells that have more hints than the other three
			unsigned short result = 0x0;
			leastCommon = h[0] & h[1] & h[2] & h[3];
			for (int i = 0; i < 4; i++)
			{
				if (h[i] & ~leastCommon)
					result |= (1 << i);
			}
			return result;


			/*if (h[0] == h[1])
			{
				if (h[2] == h[3])
				{
					if (h[1] == h[2])
					{

					}
					else
					{

					}
				}
				else
				{
					if (h[0] == h[2])
					{
						return 3;
					}
					else if (h[0] == h[3])
					{
						return 2;
					}
				}
			}
			else
			{
				if (h[2] == h[3])
				{
					if (h[0] == h[2])
					{
						return 1;
					}
					else if (h[1] == h[2])
					{
						return 0;
					}
				}
				else
				{
					
				}
			}
			return -1;*/
		}
		void find_ur(Sudoku* s, std::vector<urResult>& results) const
		{
			short ur[4];
			for (int i = 0; i< s->nn(); i++)
			{	
				//i = 1;
				auto indexerR1 = s->getIndexer(Sudoku::Indexer::row, i);
				for (int ii = i+1; ii < s->nn(); ii++)
				{		
					//ii = 5;
					auto indexerR2 = s->getIndexer(Sudoku::Indexer::row, ii);
					for (int j = 0; j < s->nn(); j++)
					{		
						//j = 1;
						ur[0] = (*indexerR1)(j);
						ur[1] = (*indexerR2)(j);
						for (int jj = j + 1; jj < s->nn(); jj++)
						{
							//jj = 2;
							ur[2] = (*indexerR1)(jj);
							ur[3] = (*indexerR2)(jj);
							if (bitcount((1 << s->getBox(ur[0])) |
								(1 << s->getBox(ur[1])) |
								(1 << s->getBox(ur[2])) |
								(1 << s->getBox(ur[3]))) == 2)
							{
								unsigned short hints[4], urhints;
								for(int k = 0;k < 4; k++)
									hints[k] = s->cells[ur[k]].hints;

								unsigned short addis = additionals(hints, urhints);


								if(bitcount(urhints)==2)
									switch (addis)
									{
									case 1:
									case 2:
									case 4:
									case 8:
										// type 1
										results.push_back({ (1u << i) | (1u << ii), (1u << j) | (1u << jj), urhints, 1, {ur[get_last_set(addis)]}, urhints });
										break;
									case 3:
									case 5:
									case 10:
									case 12:
										// type 2
										if (bitcount(bitunion(bitindex(hints, addis, 4)) & ~urhints) == 1)
										{
											unsigned short commonAddis = bitintersection(bitindex(hints, addis, 4)) & ~urhints;
											if (commonAddis)
											{
												std::set<short> visiblePositions;
												visible(s, bitindex(ur, addis, 4), visiblePositions, false);
												results.push_back({ (1u << i) | (1u << ii), (1u << j) | (1u << jj), urhints, 2, visiblePositions, commonAddis });
											}
										}
										break;
									}
							}
						}
					}
				}
			}
		}
		void handle_ur(Sudoku* s, urResult urr, bool& changed) const
		{
			// highlight the board
			for (int i = 0; i < s->nn(); i++)
				if (urr.urrows & (1u << i))
				{
					s->markers[{-1, i, 0}] = { 1 };
					for (int j = 0; j < s->nn(); j++)
						if (urr.urclms & (1u << j))
						{
							s->markers[{j, -1, 0}] = { 1 };
							for (int h = 0; h < s->nn(); h++)
								if (urr.urhints & (1u << h))
								{
									s->markers[{j, i, h + 1}] = { 0 };
								}
						}
				}

			// erase
			int nelims = 0;
			for (auto c : urr.ercells)
			{
				uint16_t& cellhints = s->cells[c].hints;
				if (!cellhints)
					continue;
				nelims += bitcount(cellhints);
				cellhints &= ~urr.erhints;
				nelims -= bitcount(cellhints);
				for (int h = 0; h < s->nn(); h++)
					if (urr.erhints & (1u << h))
						s->markers[{s->getClm(c), s->getRow(c), h + 1}] = { 2 };
			}
			std::cout << " -> eliminate " << nelims << " occurences";
			for (auto c : urr.ercells)
			{
				// set
				uint16_t& cellhints = s->cells[c].hints; 
				if (!cellhints)
					continue;
				if ((cellhints & (cellhints - 1)) == 0) /// only one bit set?
				{
					std::cout << " -> set " << getTextFromBits(cellhints, 1) << " at r" << s->getRow(c) + 1 << "c" << s->getClm(c) + 1;
					s->markers[{s->getClm(c), s->getRow(c), 0}] = { 0 };

					s->cells[c].value = cellhints;
					eliminateHints(s, cellhints,
						s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(c)),
						s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(c)),
						s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(c)));
					changed = true;
				}
			}

		}

	public: 
		void apply(Sudoku* s, int maxorder) const override
		{
			s->markers.clear();
			std::vector<urResult> results;
			find_ur(s, results);
			bool hasChanged = false;
			for (auto& a : results)
			{
				std::cout << "found Unique Rectangle Type " << a.urtype << " of candidates " << getTextFromBits(a.urhints, 1) << " in " 
					<< descr(Sudoku::Indexer::row) << getTextFromBits(a.urrows, 1, "")
					<< descr(Sudoku::Indexer::clm) << getTextFromBits(a.urclms, 1, "");

				handle_ur(s, a, hasChanged);
				std::cout << std::endl;

				//if (hasChanged)
					break;
			}
		}
	};

	class BruteforceSolver : public SudokuSolver
	{
		mutable int iters = 0;
		bool brut(Sudoku* s, const std::vector<Cell>& cells, int startcell) const
		{
			iters++;
			if (startcell == s->nn() * s->nn())
			{
				// all cells successfully filled
				s->cells = cells;
				return true;
			}
			if (cells[startcell].value) return brut(s, cells, startcell + 1);
			for (int h = 0; h < s->nn(); h++)
			{
				if (cells[startcell].hints & (1 << h))
				{
					auto cells2 = cells;
					cells2[startcell].value = (1 << h);
					eliminateHints(cells2, s->nn(), (1 << h),
						s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(startcell)),
						s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(startcell)),
						s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(startcell)));

					if (brut(s, cells2, startcell + 1))
						return true;
				}
			}
			return false;
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			if (brut(s, s->cells, 0))
				std::cout << setColor(lgreen) << "Solution found! Iterations: " << iters << setColor() << std::endl;
			else
				std::cout << setColor(lred) << "There is no solution! Iterations: " << iters << setColor() << std::endl;
		}
	};
}
