#pragma once
#include <vector>
#include <deque>
#include <tuple>
#include <map>

short count_bit_o(unsigned int i)
{
	// Java: use >>> instead of >>
	// C or C++: use uint32_t
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
unsigned short count_bit(unsigned short x)
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

	protected:
		int j_;
		int n;
		int nn;
		Indexer& operator[](int j)
		{
			j_ = j;
			return *this;
		}
	};

	std::vector<Cell> cells; // rowwise, as always...
	Indexer::Type types [3] = { Indexer::Type::row, Indexer::Type::clm, Indexer::Type::box };

	struct Marker
	{
		int hint; // 0 for the cell value
		int markGroup;
	};
	std::map<std::pair<short,short>, Marker> markers;

	Sudoku(int n = 0)
	{
		init(n);
	}
	void init(int n)
	{
		n_ = n;
		nn_ = n * n;

		row_indexer.n = clm_indexer.n = box_indexer.n = n_;
		row_indexer.nn = clm_indexer.nn = box_indexer.nn = nn_;

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
	void set(int x, int y, int value)
	{
		cells[y * nn_ + x].value = 1 << (value-1);
		cells[y * nn_ + x].hints = 0;
	}
	int get(int x, int y)
	{
		unsigned short bits = cells[y * nn_ + x].value;
		return get_last_set(bits) + 1;
	}
	bool hint(int x, int y, int hint)
	{
		return (bool) (cells[y * nn_ + x].hints & (1 << (hint - 1)));
	}
	void out()
	{
		using namespace std;
		std::vector<color> cols = { turque,orange,lblue,lpurple,green,red,yellow };
		int maxiter = 0;
		for (int y = 0; y < nn_+1; y++)
		{
			cout << "   ";
			if (y % n_ == 0)
				cout << string(nn_*2+n_+(n_+1), '-') << endl << "   ";
			if (y < nn_)
			{
				for (int x = 0; x < nn_; x++)
				{
					if (x%n_ == 0)
						cout << "| ";
					int value = get(x, y);
					bool mark = markers.contains({ x,y }) && markers[{x, y}].hint == 0;
					int iter = markers[{x, y}].markGroup;
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
		std::vector<color> foreColors = { green, turque, lred };
		std::vector<color> backColors = { black, gray, blue};
		int rowhl = 0, colhl = 0;

		using namespace std;
		for (int y = 0; y < nn_+1; y++)
		{
			if (y % n_ == 0)
				cout << string(nn_ * (n_ + 2) + n_*2 + (n_ + 1), '-') << endl;
			rowhl = markers.contains({ -1,y }) ? 1 : 0;
			if (y < nn_)
			{
				for (int i = 0; i < n_; i++)
				{
					for (int x = 0; x < nn_; x++)
					{
						colhl = markers.contains({ x,-1 }) ? 1 : 0;
						if(x == 0)
							std::cout << setColor() << "|"<< setColor(white, backColors[rowhl]) <<"  ";
						else if (x % n_ == 0)
							std::cout << setColor(white, backColors[rowhl]) << "|  ";

						for (int ii = 0; ii < n_; ii++)
						{
							short h = i * n_ + ii + 1;
							bool mark = markers.contains({ x,y }) && markers[{x, y}].hint == h;
							auto& m = markers[{x, y}];
							std::cout << setColor(mark ? foreColors[m.markGroup] : white, backColors[rowhl+colhl]) <<
								(hint(x, y, h)||mark ? to_string(h) : " ") << setColor(white, backColors[rowhl]);
						}

						std::cout << "  ";
					}
					cout << setColor() << "|"  << endl;
				}
				
				if((y+1)%n_!=0)
					cout << endl;
			}
		}

		cout << endl;
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
		catch (std::exception& e)
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
			ret += count_bit(a.hints);
		return ret;
	}


private:
	struct : public Indexer {
		int operator()(int i) override
		{
			return i + j_*nn;
		}
		Type type() override { return Indexer::Type::row; }
	}row_indexer;

	struct : public Indexer {
		int operator()(int i) override
		{
			return i*nn + j_;
		}
		Type type() override { return Indexer::Type::clm; }
	}clm_indexer;

	struct : public Indexer {
		int operator()(int i) override
		{
			int base = (j_%n)*n + (j_ / n)*n*nn; // base aus box-index j
			int box = base + i%n + (i / n)*nn;
			return box;
		}
		Type type() override { return Indexer::Type::box; }
	}box_indexer;
};


int getRow(int cell, int n)
{
	return cell / (n*n); // get row index from cell index
}
int getClm(int cell, int n)
{
	return cell % (n*n); // get column index from cell index
}
int getBox(int cell, int n)
{
	return (cell / n) % n + n*(cell / (n*n*n)); // get box index from cell index
}

std::string getTextFromBits(uint16_t bits, int baseIndex)
{
	std::string ret;
	for (int value = 0; value < 16; value++)
	{
		if ((bits & (1 << value)) != 0)
			ret += (ret == "" ? "" : ", ") + std::to_string(value+ baseIndex);
	}
	return ret;
}

bool eliminateHints(Sudoku* s, uint16_t hints, Sudoku::Indexer* indexer0, 
	Sudoku::Indexer* indexer1 = nullptr, Sudoku::Indexer* indexer2 = nullptr, 
	uint16_t subindex_mask = 0xffff)
{
	bool changed_something = false;
	for (int i = 0; i < s->nn(); i++)
	{
		if (((1 << i) & subindex_mask) == 0)
			continue;

		// eliminate hints in first indexer (row propably)
		if (indexer0)
		{
			uint16_t& cellhints = s->cells[(*indexer0)(i)].hints;
			if ((cellhints & hints) != 0)
				changed_something = true;
			cellhints &= ~hints;
		}
			 

		// eliminate hints in second indexer (column propably)
		if (indexer1)
		{
			uint16_t& cellhints = s->cells[(*indexer1)(i)].hints;
			if ((cellhints & hints) != 0)
				changed_something = true;
			cellhints &= ~hints;
		}

		// eliminate hints in third indexer (box propably)
		if (indexer2)
		{
			uint16_t& cellhints = s->cells[(*indexer2)(i)].hints;
			if ((cellhints & hints) != 0)
				changed_something = true;
			cellhints &= ~hints;
		}
	}
	return changed_something;
}
void initialCalcHints(struct Sudoku* s)
{
	for (int cell = 0; cell < s->nn()*s->nn(); cell++)
	{
		if (s->cells[cell].value != 0)
		{
			eliminateHints(s, s->cells[cell].value,
				s->getIndexer(Sudoku::Indexer::Type::row, getRow(cell, s->n())),
				s->getIndexer(Sudoku::Indexer::Type::clm, getClm(cell, s->n())),
				s->getIndexer(Sudoku::Indexer::Type::box, getBox(cell, s->n())));
		}
	}
}

//--------------------------------------------------------------------------------------------
void find_naked(Sudoku* s, Sudoku::Indexer* indexer, 
	int z, int m, std::vector<uint16_t>& results, 
	uint16_t hints = 0x0, uint16_t positions = 0x0, int start_subindex = 0)
{
	if (z == 0)
	{
		results.push_back(hints);
		results.push_back(positions);
		return;
	}

	for (short i = start_subindex; i < s->nn(); i++)
	{
		if (s->cells[(*indexer)(i)].hints == 0x0)
			continue;		

		// merge hints in structure of given indexer
		uint16_t newhints = hints | s->cells[(*indexer)(i)].hints;
		if (count_bit(newhints) > m)
			continue;

		uint16_t newpositions = positions | (1 << i);

		find_naked(s, indexer, z - 1, m, results, newhints, newpositions, i+1);
	}
}
void handle_naked(Sudoku* s, Sudoku::Indexer* indexer, 
	uint16_t positions, uint16_t hints, bool& changed, int iter)
{
	// erase those hints from elsewhere in the structure
	if ((positions & (positions - 1)) == 0) /// only one bit set?
	{
		int cell = (*indexer)(get_last_set(positions));
		s->cells[cell].value = hints;
		eliminateHints(s, hints,
			s->getIndexer(Sudoku::Indexer::Type::row, getRow(cell, s->n())),
			s->getIndexer(Sudoku::Indexer::Type::clm, getClm(cell, s->n())),
			s->getIndexer(Sudoku::Indexer::Type::box, getBox(cell, s->n())));
		changed = true;
		std::cout << " -> set " << getTextFromBits(hints, 1);
		s->markers[{getClm(cell, s->n()), getRow(cell, s->n())}] = { 0,iter-1 };
	}
	else
		changed |= eliminateHints(s, hints, indexer, nullptr, nullptr, ~positions);
}

void naked(Sudoku* s, int maxorder = 8)
{
	/// are there m cells in a structure that share z <= m hints? Then remove these hints from elsewhere in the structure (here: m = z)
	s->markers.clear(); 
	std::vector<uint16_t> results; // first the hints, second sub-index positions
	bool changed = true;
	int iter = 1;
	while (changed)
	{
		changed = false;
		for (int z = 1; z <= 8; z++)
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
						std::cout << "found naked " << z << " in " << (t == Sudoku::Indexer::Type::row ? "row" : t == Sudoku::Indexer::Type::clm ? "column" : "box")
							<< j << ": hints " << getTextFromBits(results[r],1) << " on pos " << getTextFromBits(results[r + 1],0);
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

//--------------------------------------------------------------------------------------------
uint16_t getPositions(Sudoku* s, Sudoku::Indexer* indexer, uint16_t hints)
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

void find_hidden(Sudoku* s, Sudoku::Indexer* indexer, int z, int m, std::vector<uint16_t>& results, uint16_t hints = 0x0, uint16_t positions = 0x0, int start_hint = 0)
{
	if (z == 0)
	{
		results.push_back(hints);
		results.push_back(positions);
		return;
	}

	for (short h = start_hint; h < s->nn(); h++)
	{
		uint16_t thisPositions = getPositions(s, indexer, (1 << h));
		if (thisPositions == 0x0)
			continue;

		// merge positions in structure of given indexer
		uint16_t newpositions = positions | thisPositions;
		if (count_bit(newpositions) > m)
			continue;

		uint16_t newhints = hints | (1 << h);

		find_hidden(s, indexer, z - 1, m, results, newhints, newpositions, h + 1);
	}
}
void handle_hidden(Sudoku* s, Sudoku::Indexer* indexer, uint16_t positions, uint16_t hints, bool& changed, int iter)
{
	// erase every other hints from these positions
	if ((positions & (positions - 1)) == 0) /// only one bit set?
	{
		int cell = (*indexer)(get_last_set(positions));
		s->cells[cell].value = hints;
		s->cells[cell].hints = 0x0;
		eliminateHints(s, hints,
			s->getIndexer(Sudoku::Indexer::Type::row, getRow(cell, s->n())),
			s->getIndexer(Sudoku::Indexer::Type::clm, getClm(cell, s->n())),
			s->getIndexer(Sudoku::Indexer::Type::box, getBox(cell, s->n())));
		changed = true;
		std::cout << " -> set " << getTextFromBits(hints, 1);
		s->markers[{getClm(cell, s->n()), getRow(cell, s->n())}] = { 0,iter- 1 };
	}
	else
		changed |= eliminateHints(s, ~hints, indexer, nullptr, nullptr, positions);
}

void hidden(Sudoku* s, int maxorder = 8)
{
	/// are there z hints in a structure that distribute over m <= z cells?
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
						std::cout << "found hidden " << z << " in " << (t == Sudoku::Indexer::Type::row ? "row" : t == Sudoku::Indexer::Type::clm ? "column" : "box")
							<< j << ": hints " << getTextFromBits(results[r], 1) << " on pos " << getTextFromBits(results[r + 1], 0);
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

//--------------------------------------------------------------------------------------------
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
bool subtract(std::vector<bool> positions0, std::vector<bool> positions1, std::vector<bool>& out) // returns if 1 is entirely in 0
{
	for (int cell = 0; cell < positions0.size(); cell++)
	{
		if (positions0[cell] && !positions1[cell])
		{
			out[cell] = true;
		}
		if (!positions0[cell] && positions1[cell])
		{
			return false;
		}
	}
	return true;
}
int count(std::vector<bool> const& positions)
{
	int count{ 0 };
	for (const auto& a : positions)
		if (a) count++;

	return count;
}

std::vector<bool> basePositions;
void find_fish(Sudoku* s, uint16_t hints, int z, int m, std::vector<bool> positions,
	std::deque<std::pair<Sudoku::Indexer::Type, int>> sets, std::vector<std::vector<bool>>& results, 
	std::vector<std::pair<Sudoku::Indexer::Type, int>>& locations, bool cover = false, int start_j = 0)
{
	if (z == 0)
	{
		if (cover)
		{
			// do the cover sets cover all given hints in the base sets?
			auto temp = std::vector<bool>(s->nn()*s->nn(), false);
			if (subtract(positions, basePositions, temp))
			{
				if (temp.size() != 0)
				{
					// mark the hints
					s->markers.clear();
					for (int i = 0; i < positions.size(); i++)
						if (positions[i])
							s->markers[{i% s->nn(), i / s->nn()}] = { get_last_set(hints) + 1,0 };
					for (int i = 0; i < basePositions.size(); i++)
						if (basePositions[i])
							s->markers[{i% s->nn(), i / s->nn()}] = { get_last_set(hints) + 1,1 };
					for (int i = 0; i < temp.size(); i++)
						if (temp[i])
							s->markers[{i% s->nn(), i / s->nn()}] = { get_last_set(hints) + 1,2 };
					for (auto& se : sets)
						if (se.first == Sudoku::Indexer::row)
							s->markers[{-1, se.second}] = { 0,0 };
						else
							s->markers[{se.second, -1}] = { 0,0 };
				}

				results.push_back(temp);
				for (int zz = 0; zz < 2 * m; zz++)
				{
					locations.push_back(sets[zz]);
				}
			}
			return;
		}
		else
		{
			// Reset the values. The base sets are found, Now look for cover sets.
			basePositions = positions;
			//basePositions.flip();
			positions = std::vector<bool>(s->nn()*s->nn(), false);

			z = m;
			cover = true;
			start_j = 0;
		}
	}
	// for every structure
	for (Sudoku::Indexer::Type t : s->types)
	{
		if (t == Sudoku::Indexer::Type::box) // TODO fishes with boxes are not (yet) supported
			continue;
		if (cover && t == sets[0].first) // dont allow base- and cover-sets to be of same type
			continue;
		for (int j = start_j; j < s->nn(); j++)
		{
			// get indexer
			Sudoku::Indexer* indexer = s->getIndexer(t, j);

			// if there is no given hint in this structure or a hint is also in an already chosen structure, this is invalid
			uint16_t thisPositions = getPositions(s, indexer, hints);
			if (thisPositions == 0x0 || intersects(s, indexer, positions, hints))
				continue;

			// merge positions in structure of given indexer
			auto newpositions = merge(s, indexer, positions, hints);
			/*if (count(newpositions) > m)
				continue;*/

			auto newSets = sets;
			newSets.push_back({ t, j });

			// go on to find next structure
			find_fish(s, hints, z - 1, m, newpositions, newSets, results, locations, cover, j+1);

			newSets.pop_back();
		}
	}
}

void handle_fish(Sudoku* s, std::vector<bool> positions, uint16_t hints, bool& changed)
{
	// erase every other hints from these positions
	for (int cell = 0; cell < s->nn()*s->nn(); cell++)
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
				std::cout << " -> set " << getTextFromBits(cellhints, 1) << " at r" << getRow(cell, s->n()) << "c" << getClm(cell, s->n());
				s->markers[{getClm(cell, s->n()),getRow(cell, s->n())}] = { 0,0 };

				s->cells[cell].value = cellhints;
				eliminateHints(s, cellhints,
					s->getIndexer(Sudoku::Indexer::Type::row, getRow(cell, s->n())),
					s->getIndexer(Sudoku::Indexer::Type::clm, getClm(cell, s->n())),
					s->getIndexer(Sudoku::Indexer::Type::box, getBox(cell, s->n())));
				changed = true;
			}
		}
	}
}

void fish(Sudoku* s, int maxorder = 8)
{
	std::string fishnames [] = { "X-Wing", "Swordfish", "Jellyfish" };
	s->markers.clear();
	std::vector<std::vector<bool>> results; // first the hints, second sub-index positions
	std::vector<std::pair<Sudoku::Indexer::Type, int>> locs;

	bool changed = true;
	while (changed)
	{
		changed = false;
		for (int z = 2; z <= maxorder; z++)
		{
			if (z > 4)
			{
				std::cout << "Maximal fish order is 4. For 9x9 Sudoku this covers all cases. See http://hodoku.sourceforge.net/de/tech_fishb.php#bf5" << std::endl;
				break;
			}
			for (short h = 0; h < s->nn(); h++)
			{
				find_fish(s, (1 << h), z, z, std::vector<bool>(s->nn()*s->nn(), false), {}, results, locs);
				for (int r = 0; r < results.size() ; r += 1)
				{
					std::cout << "found " << fishnames[z-2] << " ("<<z << ") of hint " << (h + 1) << " in ";
					for (int zz = 0; zz < 2*z; zz++)
					{
						if (zz == z)
							std::cout << " ";
						auto t = locs[r*2*z+zz].first;
						std::cout << (t == Sudoku::Indexer::Type::row ? "r" : t == Sudoku::Indexer::Type::clm ? "c" : "b") << locs[r * 2 * z + zz].second;
					}
					std::cout << " delete: " << count(results[r]);
					handle_fish(s, results[r], (1 << h), changed);
					std::cout << std::endl;
					if(changed)
						return;
				}
				results.clear();
				locs.clear();
			}
			
		}
		break;
	}
}

