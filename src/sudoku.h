#pragma once
#include "dynamic_bitset.h"

#include <vector>
#include <deque>
#include <tuple>
#include <map>
#include <set>
#include <array>
#include <assert.h>
#include <sstream>
#include <limits>
#include <algorithm>

namespace sudoku
{
	using index_t = unsigned short; // Must be able to address all n^2 cells.
	using candidates_t = uint16_t;  // Candidates are the single bits.
	using positions_t = uint16_t;   // Positions [0,n) are the single bits.

	// Returns the number of bits in x.
	constexpr uint16_t bitcount(uint16_t x)
	{
		// TODO use std::popcount() once it's available
		x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
		x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
		x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
		return x;
	}

	// Returns the index of the last set bit (eg.: 4->2, 6->1).
	constexpr int get_last_set(uint16_t bits)
	{
		if (bits == 0)
			return -1;
		int value = 0;
		while ((bits & 1) == 0){
			bits >>= 1;
			value++;
		}
		return value;
	}


	// Class for a single Sudoku cell.
	struct Cell
	{
		candidates_t value;       // The value is the index of the single bit set
		candidates_t candidates;  // Zero or many set bits possible
	};

	// Class for a full Sudoku.
	class Sudoku
	{
		index_t n_; // The block size, usually 3
		index_t nn_; // The square of the block size = sudoku size
		inline void init(index_t n)
		{
			n_ = n;
			nn_ = n * n;

			for (index_t j = 0; j < nn_; j++)
			{
				row_indexer.push_back({ j, n_, nn_ });
				clm_indexer.push_back({ j, n_, nn_ });
				box_indexer.push_back({ j, n_, nn_ });
			}

			cells.resize(nn_ * nn_);
			for (index_t i = 0; i < nn_ * nn_; i++)
			{
				cells[i].value = 0;
				cells[i].candidates = 0xffff >> (16 - nn_);
			}
		}

		// A helper class to print the cell values of a Sukoku.
		struct ValueRepresentation{
		private:
			friend class Sudoku;
			Sudoku* s;
			ValueRepresentation(Sudoku* s):s{s}{}
			std::ostream& print(std::ostream& os) const {
				std::vector<color> cols ={turque,orange,lblue,lpurple,green,red,yellow,purple,blue,gray,gray};
				index_t maxiter = 0;
				for(index_t y = 0; y < s->nn_ + 1; y++)
				{
					os << "   ";
					if(y % s->n_ == 0)
						os << std::string(s->nn_ * 2 + s->n_ + (s->n_ + 1),'-') << std::endl << "   ";
					if(y < s->nn_)
					{
						for(index_t x = 0; x < s->nn_; x++)
						{
							if(x % s->n_ == 0)
								os << "| ";
							int value = s->getCellValue(x,y);
							bool mark = s->markers.contains({x,y,0});
							index_t iter = s->markers[{x,y,0}].markGroup;
							maxiter = std::max(maxiter,iter);
							if(value != 0)
								os << setColor(mark ? cols[iter] : white) << value << setColor(white) << " ";
							else
								os << "  ";
						}
						os << "|" << std::endl;
					}
				}
				os << std::endl << "Iterations: ";
				for(index_t i = 0; i <= maxiter; i++)
				{
					os << setColor(cols[i]) << i << setColor() << (i == maxiter ? "." : ", ");
				}
				return os << std::endl;
			}
		public:
			friend std::ostream& operator<< (std::ostream& os, ValueRepresentation const& cr){
				return cr.print(os);
			}
		};

		// A helper class to print the candidates of a Sukoku.
		struct CandidateRepresentation{
		private:
			friend class Sudoku;
			Sudoku* s;
			CandidateRepresentation(Sudoku* s):s{s}{}
			std::ostream& print(std::ostream& os) const {
				std::vector<color> foreColors ={turque, orange, red, lgreen};
				std::vector<color> backColors ={black, gray, blue, purple};
				int rowhl = 0,colhl = 0,boxhl = 0;
				bool marked = false;
				for(int y = 0; y < s->nn_ + 1; y++)
				{
					if(y % s->n_ == 0)
						os << std::string(s->nn_ * (s->n_ + 2) + s->n_ * 2 + (s->n_ + 1),'-') << std::endl;
					rowhl = s->markers.contains({-1,y,0}) ? 1 : 0;
					if(y < s->nn_)
					{
						for(int iy = 0; iy < s->n_; iy++)
						{
							for(int x = 0; x < s->nn_; x++)
							{
								colhl = s->markers.contains({x,-1,0}) ? 1 : 0;
								boxhl = s->markers.contains({-1,-1,s->getBox(s->nn_*y+x)}) ? 1 : 0;
								if(x == 0)
									os << setColor() << "|" << setColor(white,backColors[rowhl]) << "  ";
								else if(x % s->n_ == 0)
									os << setColor(white,backColors[rowhl]) << "|  ";

								for(int ix = 0; ix < s->n_; ix++)
								{
									short h = iy * s->n_ + ix + 1;
									bool mark = s->markers.contains({x,y,h});
									marked |= mark;
									auto& m = s->markers[{x,y,h}];
									bool forced = !s->hasCandidate(x,y,h) && mark;
									os << setColor(mark ? foreColors[m.markGroup] : white,forced ? yellow : backColors[rowhl + colhl + boxhl]) <<
										((s->hasCandidate(x,y,h) || mark) ? std::to_string(h) : " ") << setColor(white,backColors[rowhl]);
								}

								os << "  ";
							}
							os << setColor() << "|" << std::endl;
						}

						if((y + 1) % s->n_ != 0)
							os << std::endl;
					}
				}
				if(marked)
					os << "(" << setColor(foreColors[0]) << "base" << setColor() << ", " << setColor(foreColors[1]) << "cover" << setColor() << ", " << setColor(foreColors[2]) << "eliminated" << setColor() << ", " << setColor(foreColors[3]) << "fin" << setColor() << ")" << std::endl;
				return os;
			}
		public:
			friend std::ostream& operator<< (std::ostream& os,CandidateRepresentation const& cr){
				return cr.print(os);
			}
		};

	public:
		// An indexer allows to access the cells within a substructure (a house) of the Sudoku.
		// For the standard Sudokus three kinds of houses are important: rows, columns and boxes.
		struct Indexer
		{
			enum class Type : size_t
			{
				row=0,
				clm,
				box
			};
			virtual index_t operator ()(index_t i) const = 0;
			virtual Type type() const = 0;
			index_t j() const
			{
				return j_;
			}
			Indexer(index_t j, index_t n, index_t nn) : j_{ j }, n{ n }, nn{ nn }{}
		protected:
			index_t j_;
			index_t n;
			index_t nn;
		};

		std::vector<Cell> cells; // rowwise, as always...
		inline static const Indexer::Type types[3] = { Indexer::Type::row, Indexer::Type::clm, Indexer::Type::box };

		struct Marker
		{
			int markGroup;
		};
		std::map<std::tuple<index_t, index_t, candidates_t>, Marker> markers; // x,y,hint(0=value)

		ValueRepresentation valueRep;
		CandidateRepresentation candiRep;

		Sudoku(int n = 0):valueRep(this),candiRep(this)
		{
			init(n);
		}
		Sudoku(const Sudoku& s) : Sudoku(s.n_){
			cells=s.cells;
		}
		Sudoku& operator=(const Sudoku& s) = delete;

		int n() const {
			return n_;
		}
		int nn() const {
			return nn_;
		}
		int nnn() const {
			return nn_*n_;
		}

		// Sets the value at cell coordinate [x,y].
		void setCellValue(index_t x, index_t y, index_t value)
		{
			cells[y * nn_ + x].value = 1ul << (value - 1);
			cells[y * nn_ + x].candidates = 0;
		}

		// Returns the value at cell coordinate [x,y].
		int getCellValue(index_t x, index_t y) const
		{
			return get_last_set(cells[y * nn_ + x].value) + 1;
		}

		// Returns if cell [x,y] allows candidate candi.
		bool hasCandidate(index_t x, index_t y, index_t candi) const
		{
			return (bool)(cells[y * nn_ + x].candidates & (1ul << (candi - 1)));
		}


		// Loads from a text file. Returns if an error occured.
		bool loadFromFile(std::string filename)
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
							setCellValue(x, y, value);
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


		// Returns if all cells of this sudoku have values.
		bool isSolved() const
		{
			for (auto& a : cells)
				if (a.value == 0)
					return false;
			return true;
		}

		// Returns the total number of cells with value.
		int getNumValues() const
		{
			int ret = 0;
			for (auto& a : cells)
				ret += (a.value != 0);
			return ret;
		}

		// Returns the total number of open candidates in all cells.
		int getNumCandidates() const
		{
			int ret = 0;
			for (auto& a : cells)
				ret += bitcount(a.candidates);
			return ret;
		}

		// Returns a pointer to the j'th indexer of type t.
		const Indexer* getIndexer(Indexer::Type t, index_t j) const
		{
			if (t == Indexer::Type::row)
				return &row_indexer[j];
			else if (t == Indexer::Type::clm)
				return &clm_indexer[j];
			else if (t == Indexer::Type::box)
				return &box_indexer[j];
			else return nullptr;
		}

		// Returns the row index from cell index.
		index_t getRow(index_t cell) const
		{
			return cell / (nn_); 
		}
		// Returns the column index from cell index.
		index_t getClm(index_t cell) const
		{
			return cell % (nn_);
		}
		// Returns the box index from cell index.
		index_t getBox(index_t cell) const
		{
			return (cell / n_) % n_ + n_ * (cell / (nn_*n_));
		}
		// Returns the index of the given house type from cell index.
		index_t getHouse(Indexer::Type t, index_t cell) const
		{
			switch(t)
			{
			case Indexer::Type::row:
				return getRow(cell);
			case Indexer::Type::clm:
				return getClm(cell);
			case Indexer::Type::box:
				return getBox(cell);
			}
		}

	private:
		struct rowIndexer : public Indexer {
			using Indexer::Indexer;
			index_t operator()(index_t i) const override
			{
				return i + j_ * nn;
			}
			Type type() const override { return Indexer::Type::row; }
		};
		std::vector<rowIndexer> row_indexer;

		struct clmIndexer : public Indexer {
			using Indexer::Indexer;
			index_t operator()(index_t i) const override
			{
				return i * nn + j_;
			}
			Type type() const override { return Indexer::Type::clm; }
		};
		std::vector<clmIndexer> clm_indexer;

		struct boxIndexer: public Indexer {
			using Indexer::Indexer;
			index_t operator()(index_t i) const override
			{
				int base = (j_ % n) * n + (j_ / n) * n * nn; // base aus box-index j
				int box = base + i % n + (i / n) * nn;
				return box;
			}
			Type type() const override { return Indexer::Type::box; }
		};
		std::vector<boxIndexer> box_indexer;
	};


	/// ===================== Now come some handy algorithms that will be useful for solving. =====================

	std::string descr(Sudoku::Indexer::Type t)
	{
		return t == Sudoku::Indexer::Type::row ? "r" : t == Sudoku::Indexer::Type::clm ? "c" : "b";
	}

	// Returns a readable enumeration of the bits set.
	std::string bitsToString(uint16_t bits, int baseIndex, std::string seperator = ", ")
	{
		std::string ret;
		for (int value = 0; value < 16; value++)
		{
			if ((bits & (1ul << value)) != 0)
				ret += (ret == "" ? "" : seperator) + std::to_string(value + baseIndex);
		}
		return ret;
	}

	// Eliminate all given candidates in all given positions in all given indexers.
	int eliminateCandidates(Sudoku* s, candidates_t candidates, std::vector<const Sudoku::Indexer*> indexers, uint16_t positions = -1)
	{
		int numChanges = 0;
		for (int i = 0; i < s->nn(); i++)
		{
			if (((1ul << i) & positions) == 0)
				continue;

			for(auto& indexer : indexers)
			{
				uint16_t& cellcandidates = s->cells[(*indexer)(i)].candidates;
				if ((cellcandidates & candidates) != 0)
					numChanges++;
				cellcandidates &= ~candidates;
			}
		}
		return numChanges;
	}

	// Eliminates all candidates in the same house as a value.
	void initialCalcCandidates(Sudoku* s)
	{
		for (index_t cell = 0; cell < s->nn() * s->nn(); cell++)
		{
			if (s->cells[cell].value != 0)
			{
				eliminateCandidates(s,s->cells[cell].value,
					{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
					 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
					 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell))});
			}
		}
	}

	// Returns if any of the given candidates are in the intersection of the two given houses.
	bool intersects(Sudoku* s, const Sudoku::Indexer* indexer0, const Sudoku::Indexer* indexer1, candidates_t candidates)
	{
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			int index0 = (*indexer0)(i0);
			for (int i1 = 0; i1 < s->nn(); i1++)
			{
				if (index0 == (*indexer1)(i1))
				{
					if ((s->cells[index0].candidates & candidates) != 0x0)
						return true;
				}
			}
		}
		return false;
	}

	// Returns if any of the given candidates are in the intersections of the house and the positions.
	bool intersects(Sudoku* s, const Sudoku::Indexer* indexer0, dyn_bs const& positions, candidates_t candidates)
	{
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			int index0 = (*indexer0)(i0);
			if (positions[index0])
			{
				if ((s->cells[index0].candidates & candidates) != 0x0)
					return true;
			}
		}
		return false;
	}

	// check if any of the given candidates are in intersections between structures
	dyn_bs merge(Sudoku* s, const Sudoku::Indexer* indexer0, dyn_bs positions, candidates_t candidates)
	{
		dyn_bs ret = positions;
		for (int i0 = 0; i0 < s->nn(); i0++)
		{
			auto index0 = (*indexer0)(i0);
			if ((s->cells[index0].candidates & candidates) != 0x0)
			{
				ret.set(index0, true);
			}
		}
		return ret;
	}

	// Unsets all positions that are not visible from any of the given cells.
	void visible(const Sudoku* s, dyn_bs const& cells, dyn_bs& out_positions)
	{
		for (int j = 0; j < s->nn() * s->nn(); j++)
		{
			if (!cells[j])
				continue;
			int cr = s->getRow(j);
			int cc = s->getClm(j);
			int cb = s->getBox(j);
			for (int i = 0; i < s->nn() * s->nn(); i++)
			{
				out_positions.and_assign(i, ((cr == s->getRow(i)) || (cc == s->getClm(i)) || (cb == s->getBox(i))));
			}
		}
	}

	// Sets positions to the cells that are visible by all given cells.
	void visible(const Sudoku* s, std::vector<index_t> const& cells, std::vector<index_t>& out_positions, bool includeCells)
	{
		assert(cells.size() <= 16);
		std::map<index_t, index_t> os;
		for (int j = 0; auto c : cells)
		{
			for (auto t : s->types)
			{
				auto indexer = s->getIndexer(t, s->getHouse(t, c));
				for (int i = 0; i < s->nn(); i++)
					os[(*indexer)(i)] |= (1ul << j);
			}
			j++;
		}
		for (auto& a : os)
			if (bitcount(a.second) == cells.size() && (std::find(std::begin(cells),std::end(cells),a.first)==std::end(cells) || includeCells))
				out_positions.push_back(a.first);
	}

	// Returns if set 1 is entirely in set 0.
	bool subtract(dyn_bs const& positions0, dyn_bs const& positions1, dyn_bs& only0, dyn_bs& only1, bool breakOnOnly1)
	{
		bool contained = true;
		for(int cell = 0; cell < positions0.size(); cell++)
		{
			if(positions0[cell] && !positions1[cell])
			{
				only0.set(cell);
			}
			else if(!positions0[cell] && positions1[cell])
			{
				only1.set(cell);
				if(breakOnOnly1)
					return false;
				contained = false;
			}
		}
		return contained;
	}

	// sets excactly those positions0 that are not part of indexer1 but have at least one given hint.
	void subtract(Sudoku* s, const Sudoku::Indexer* indexer0, const Sudoku::Indexer* indexer1, uint16_t& positions0, uint16_t hintmask = 0xffff)
	{
		positions0 = 0x0;
		for (int i = 0; i < s->nn(); i++)
		{
			auto cell = (*indexer0)(i);
			if (s->cells[cell].candidates & hintmask)
			{
				uint16_t p = (1ul << i);
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

	// Returns the subindices of those cells, that contain all given candidates.
	positions_t containingAll(Sudoku* s, const Sudoku::Indexer* indexer, candidates_t candidates)
	{
		positions_t positions = 0;
		for (int i = 0; i < s->nn(); i++)
		{
			uint16_t temp = s->cells[(*indexer)(i)].candidates;
			if ((temp & candidates) == candidates)
				positions |= (1ul << i);
		}
		return positions;
	}

	// get the subindices of those cells, that contain any of the given candidates.
	uint16_t containingAny(Sudoku* s, const Sudoku::Indexer* indexer, candidates_t candidates)
	{
		uint16_t positions = 0x0;
		for (int i = 0; i < s->nn(); i++)
		{
			uint16_t temp = s->cells[(*indexer)(i)].candidates;
			if (temp & candidates)
				positions |= (1ul << i);
		}
		return positions;
	}

	// returns the indexer of the second house, all the given positions of the given indexer are in, if there is one.
	const Sudoku::Indexer* secondHouse(Sudoku* s, const Sudoku::Indexer* indexer, positions_t positions)
	{
		uint16_t secondHouse1 = 0xffff,secondHouse2 = 0xffff;
		for(int i = 0; i < s->nn(); i++)
		{
			if((positions & (1ul << i)) != 0)
			{
				switch(indexer->type())
				{
				case Sudoku::Indexer::Type::row: [[fallthrough]];
				case Sudoku::Indexer::Type::clm:
					secondHouse1 &= (1ul << s->getBox((*indexer)(i)));
					secondHouse2 = 0x0;
					break;
				case Sudoku::Indexer::Type::box:
					secondHouse1 &= (1ul << s->getRow((*indexer)(i)));
					secondHouse2 &= (1ul << s->getClm((*indexer)(i)));
					break;
				}

				if(!secondHouse1 && !secondHouse2)
					return nullptr;
			}
		}
		switch(indexer->type())
		{
		case Sudoku::Indexer::Type::row:
		case Sudoku::Indexer::Type::clm:
			return s->getIndexer(Sudoku::Indexer::Type::box,get_last_set(secondHouse1));
		case Sudoku::Indexer::Type::box:
			return secondHouse1 ? s->getIndexer(Sudoku::Indexer::Type::row,get_last_set(secondHouse1)) : s->getIndexer(Sudoku::Indexer::Type::clm,get_last_set(secondHouse2));;
		}
	}

	// Returns the image of the bit positions under arr[].
	template<index_t max, typename T>
	std::vector<T> bitindex(T arr[], candidates_t candidates)
	{
		std::vector<T> ret;
		for (index_t i = 0; i < max; i++)
		{
			if (candidates & (1ul << i))
				ret.push_back(arr[i]);
		}
		return ret;
	}

	// Returns the intersection of all given bitmasks.
	candidates_t bitintersection(std::vector<candidates_t> hs)
	{
		candidates_t ret = std::numeric_limits<int>::max();
		for (const auto& s : hs)
			ret &= s;
		return ret;
	}

	// Returns the union of all given bitmasks.
	candidates_t bitunion(std::vector<candidates_t> hs)
	{
		candidates_t ret = 0x0;
		for (const auto& s : hs)
			ret |= s;
		return ret;
	}


	/// ===================== In the following the different solving techniques are implemented. =====================

	class SudokuSolver
	{
	public:
		virtual void apply(Sudoku* s, int maxorder) const = 0;
	};

	class LockedSolver : public SudokuSolver
	{
		struct lockedResult
		{
			const Sudoku::Indexer* indexer;
			uint16_t candidates;
			uint16_t positions;
		};
	protected:
		void find_locked(Sudoku* s, const Sudoku::Indexer* indexer, std::vector<lockedResult>& results) const
		{
			for (short h = 0; h < s->nn(); h++)
			{
				uint16_t thisPositions = containingAll(s, indexer, (1ul << h));
				if (thisPositions == 0x0)
					continue;

				// watch out for a second house all occurences of hint h are contained in.
				if (auto indexer2 = secondHouse(s, indexer, thisPositions))
				{
					uint16_t locked_positions;
					subtract(s, indexer2, indexer, locked_positions, 1ul << h);
					if(locked_positions)
						results.push_back({ indexer2, 1u << h, locked_positions });
				}
			}
		}
		void handle_locked(Sudoku* s, const lockedResult& r, bool& changed, int iter) const
		{
			// erase every other candidates from these positions
			if ((r.positions & (r.positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*r.indexer)(get_last_set(r.positions));
				s->cells[cell].value = r.candidates;
				s->cells[cell].candidates = 0x0;
				eliminateCandidates(s,r.candidates,
									{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell))});
				changed = true;
				std::cout << " -> set " << bitsToString(r.candidates, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateCandidates(s,r.candidates,{r.indexer},r.positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there z candidates in a structure that distribute over m <= z cells?
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
						const Sudoku::Indexer* indexer = s->getIndexer(t, j);
						find_locked(s, indexer, results);
						for (auto& r : results)
						{
							std::cout << "found Locked Candidates in " << descr(t) << j + 1ul << "/" << descr(r.indexer->type()) << r.indexer->j() + 1
								<< ": candidates " << bitsToString(r.candidates, 1) << " on pos " << bitsToString(r.positions, 1);
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
		void find_hidden(Sudoku* s, const Sudoku::Indexer* indexer, int z, int m, std::vector<uint16_t>& results,
			uint16_t candidates = 0x0, uint16_t positions = 0x0, int start_hint = 0) const
		{
			if (z == 0)
			{
				// check if there are candidates to delete:
				uint16_t delcands = 0x0;
				for (int i = 0; i < s->nn(); i++)
					if (positions & (1ul << i))
						delcands |= s->cells[(*indexer)(i)].candidates;
				delcands &= ~candidates;
				if (m == 1 || (bool)delcands)
				{
					results.push_back(candidates);
					results.push_back(positions);
				}
				return;
			}

			for (short h = start_hint; h < s->nn(); h++)
			{
				uint16_t thisPositions = containingAll(s, indexer, (1ul << h));
				if (thisPositions == 0x0)
					continue;

				// merge positions in house of given indexer
				uint16_t newpositions = positions | thisPositions;
				auto nposis = bitcount(newpositions);
		
				if (nposis > m)
					continue;

				find_hidden(s, indexer, z - 1, m, results, candidates | (1ul << h), newpositions, h + 1);
			}
		}
		void handle_hidden(Sudoku* s, const Sudoku::Indexer* indexer, uint16_t positions, uint16_t candidates, bool& changed, int iter) const
		{
			// erase every other candidates from these positions
			if ((positions & (positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*indexer)(get_last_set(positions));
				s->cells[cell].value = candidates;
				s->cells[cell].candidates = 0x0;
				eliminateCandidates(s, candidates,
									{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell))});
				changed = true;
				std::cout << " -> set " << bitsToString(candidates, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateCandidates(s,~candidates,{indexer},positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there z candidates in a structure that distribute over m <= z cells?
			std::string subsetNames[] = { "Single", "Pair", "Triple", "Quadrupel", "Quintuple", "Sextuple", "Septuple", "Octuple" };
			s->markers.clear();
			std::vector<uint16_t> results; // first the candidates, second sub-index positions
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
							const Sudoku::Indexer* indexer = s->getIndexer(t, j);
							find_hidden(s, indexer, z, z, results);
							for (int r = 0; r < results.size() / 2; r += 2)
							{
								std::cout << "found Hidden " << subsetNames[z-1] << " in " << descr(t)
									<< j + 1ul << ": candidates " << bitsToString(results[r], 1) << " on pos " << bitsToString(results[r + 1], 1);
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
		void find_naked(Sudoku* s, const Sudoku::Indexer* indexer,
			int z, int m, std::vector<uint16_t>& results,
			uint16_t candidates = 0x0, uint16_t positions = 0x0, int start_j = 0) const
		{
			if (z == 0)
			{
				// check if there are candidates to delete:
				uint16_t delposis = containingAny(s, indexer, candidates) & ~positions;
				if (m == 1 || (bool)delposis)
				{
					results.push_back(candidates);
					results.push_back(positions);
				}
				return;
			}

			for (short j = start_j; j < s->nn(); j++)
			{
				if (s->cells[(*indexer)(j)].candidates == 0x0)
					continue;

				// merge candidates in structure of given indexer
				uint16_t newcandidates = candidates | s->cells[(*indexer)(j)].candidates;
				auto ncandidates = bitcount(newcandidates);
				if (ncandidates > m)
					continue;

				find_naked(s, indexer, z - 1, m, results, newcandidates, positions | (1ul << j), j + 1);
			}
		}
		void handle_naked(Sudoku* s, const Sudoku::Indexer* indexer,
			uint16_t positions, uint16_t candidates, bool& changed, int iter) const
		{
			// erase those candidates from elsewhere in the structure
			if ((positions & (positions - 1)) == 0) /// only one bit set?
			{
				int cell = (*indexer)(get_last_set(positions));
				s->cells[cell].value = candidates;
				eliminateCandidates(s, candidates,
									{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
									 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell))});
				changed = true;
				std::cout << " -> set " << bitsToString(candidates, 1);
				s->markers[{s->getClm(cell), s->getRow(cell), 0}] = { iter - 1 };
			}
			else
			{
				int nelims = eliminateCandidates(s,candidates,{indexer},~positions);
				changed |= (bool)nelims;
				std::cout << " -> eliminate " << nelims << " occurences";
			}
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			/// are there m cells in a structure that share z <= m candidates? Then remove these candidates from elsewhere in the structure (here: m = z)
			std::string subsetNames[] = { "Single", "Pair", "Triple", "Quadrupel", "Quintuple", "Sextuple", "Septuple", "Octuple" };
			s->markers.clear();
			std::vector<uint16_t> results; // first the candidates, second sub-index positions
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
							const Sudoku::Indexer* indexer = s->getIndexer(t, j);
							find_naked(s, indexer, z, z, results);
							for (int r = 0; r < results.size() / 2; r += 2)
							{
								std::cout << "found Naked " << subsetNames[z - 1] << " in " << descr(t)
									<< j + 1ul << ": candidates " << bitsToString(results[r], 1) << " on pos " << bitsToString(results[r + 1], 1);
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
		static const bool only_homogenous_fish_bases = false;
		static const bool only_exclusive_fish_sets = false;

		struct fishResult
		{
			dyn_bs positions;
			int finns;
		};

	protected:
		mutable dyn_bs basePositions;

		void find_fish(Sudoku* s, candidates_t candidates, int z, int m, dyn_bs positions,
			std::deque<std::pair<Sudoku::Indexer::Type, int>> sets, std::vector<fishResult>& results,
			std::vector<std::pair<Sudoku::Indexer::Type, int>>& locations, bool finned, bool franken, int& ncombis,
			const int maxnum = -1, bool cover = false, std::array<int, 3> start_j ={0,0,0}, uint32_t baseSets = 0x0) const
		{
			if(z == 0)
			{
				if(cover)
				{
					ncombis++;
					// Do the cover sets cover all given candidates in the base sets?
					dyn_bs onlyCover, onlyBase;
					onlyCover.resize(s->nn() * s->nn(), false);
					onlyBase.resize(s->nn() * s->nn(), false);
					if(subtract(positions, basePositions, onlyCover, onlyBase, !finned))
					{
						// Select only the cover candidates, that are visible by all fins
						auto temp = onlyCover;
						visible(s, onlyBase, onlyCover);

						if(onlyCover.any())
						{
							visible(s, onlyBase, temp);

							// mark the candidates, override the old markers
							s->markers.clear();
							for(int i = 0; i < positions.size(); i++)
								if(positions[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(candidates) + 1}] ={1}; // cover candidates
							for(int i = 0; i < basePositions.size(); i++)
								if(basePositions[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(candidates) + 1}] ={0}; // base candidates
							for(int i = 0; i < onlyBase.size(); i++)
								if(onlyBase[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(candidates) + 1}] ={3}; // fins
							for(int i = 0; i < onlyCover.size(); i++)
								if(onlyCover[i])
									s->markers[{i% s->nn(), i / s->nn(), get_last_set(candidates) + 1}] ={2}; // deletions
							for(auto& se : sets)
								if(se.first == Sudoku::Indexer::Type::row)
									s->markers[{-1, se.second, 0}] ={0};
								else if(se.first == Sudoku::Indexer::Type::clm)
									s->markers[{se.second, -1, 0}] ={0};
								else
									s->markers[{-1, -1, se.second}] ={0};

							results.push_back({onlyCover,onlyBase.popcount()});
							for(int zz = 0; zz < 2 * m; zz++)
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
					positions.fill(false);

					z = m;
					cover = true;
					start_j[0] = start_j[1] = start_j[2] = 0;
				}
			}

			// for every structure
			for(index_t ti = -1; Sudoku::Indexer::Type t : s->types)
			{
				++ti;
				if(!franken && t == Sudoku::Indexer::Type::box) // without franken, dont consider boxes
					continue;
				if (only_homogenous_fish_bases && !cover && !sets.empty() && t != sets[0].first) // dont allow base sets to be of more than one type
					continue;
				if (only_exclusive_fish_sets && cover && t == sets[0].first) // dont allow base- and cover-sets to be of same type
					continue;
				for(index_t j = start_j[ti]; j < s->nn(); j++)
				{
					if(cover && (baseSets & (1ul << (j + ti * s->nn()))))
						continue;

					// get indexer
					const Sudoku::Indexer* indexer = s->getIndexer(t, j);

					// if there is no given hint in this structure or a hint is also in an already chosen structure, this is invalid
					positions_t thisPositions = containingAll(s, indexer, candidates);
					if(thisPositions == 0x0 || intersects(s, indexer, positions, candidates))
						continue;

					// merge positions in structure of given indexer
					auto newpositions = merge(s, indexer, positions, candidates);
					/*if (count(newpositions) > m)
						continue;*/

					auto newSets = sets;
					newSets.push_back({t, j});
					start_j[ti] = j;

					// go on to find next structure
					find_fish(s, candidates, z - 1, m, newpositions, newSets, results, locations, finned, franken, ncombis, maxnum, cover, start_j, cover ? baseSets : baseSets | (1ul << (j+(size_t)t*s->nn())));
					newSets.pop_back();

					if(results.size() == maxnum)
						return;
				}
			}
		}

		void handle_fish(Sudoku* s, dyn_bs positions, uint16_t candidates, bool& changed) const
		{
			// erase every other candidates from these positions
			std::cout << " -> eliminate " << positions.popcount() << " occurences";
			for(int cell = 0; cell < s->nn() * s->nn(); cell++)
			{
				if(positions[cell])
				{
					// erase
					uint16_t& cellcandidates = s->cells[cell].candidates;
					if((cellcandidates & ~candidates) != 0)
						changed = true;
					cellcandidates &= ~candidates;

					// set
					if((cellcandidates & (cellcandidates - 1)) == 0) // only one bit set?
					{
						std::cout << " -> set " << bitsToString(cellcandidates, 1) << " at r" << s->getRow(cell) + 1ul << "c" << s->getClm(cell) + 1;
						s->markers[{s->getClm(cell), s->getRow(cell), 0}] ={0};

						s->cells[cell].value = cellcandidates;
						eliminateCandidates(s, cellcandidates,
											{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(cell)),
											 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(cell)),
											 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(cell))});
						changed = true;
					}
				}
			}
		}


	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			const std::string fishnames[] ={"X-Wing", "Swordfish", "Jellyfish"};
			s->markers.clear();
			std::vector<fishResult> results;
			std::vector<std::pair<Sudoku::Indexer::Type, int>> locs;

			bool changed = true;
			while(changed)
			{
				changed = false;
				for(int z = 2; z <= maxorder; z++)
				{
					if(z > 4)
					{
						std::cout << "Maximal fish order is 4. For 9x9 Sudokus this covers all cases. See http://hodoku.sourceforge.net/de/tech_fishb.php#bf5" << std::endl;
						break;
					}
					int ncombis = 0;
					int searchMS = 0;
					for(short h = 0; h < s->nn(); h++)
					{
						auto startT = std::chrono::high_resolution_clock::now();
						find_fish(s, (1ul << h), z, z, dyn_bs(s->nn() * s->nn(), false), {}, results, locs, true, true, ncombis, 1);
						searchMS += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startT).count();
						for(int r = 0; r < results.size(); r += 1)
						{
							std::cout << "Order " << z << " - search (until hit) took " << searchMS << " ms (" << ncombis << " combinations)" << std::endl;
							std::string attribute = (results[r].finns == 0 ? "" : std::to_string(results[r].finns) + "-finned "), sets;
							Sudoku::Indexer::Type oldtype;
							unsigned category = 0; unsigned flags=0;
							for(int zz = 0; zz <= 2 * z; zz++)
							{
								if(zz % z == 0)
								{
									category = std::max(std::max(category, flags&4), 2*unsigned(flags % 4 == 3)); // try to distinguish franken from mutant fish
								}
								if(zz == z)
								{
									sets += "/";
									oldtype = (Sudoku::Indexer::Type)(-1);
									flags = 0;
								}
								if(zz == 2 * z)
									break;
								auto t = locs[r * 2 * z + zz].first;
								flags |= (1ul << (size_t)t);
								sets += (t != oldtype ? descr(t) : "") + std::to_string(locs[r * 2 * z + zz].second + 1);
								oldtype = t;
							}
							attribute = attribute + (category == 0 ? "" : category == 1 ? "Franken " : "Mutant ");
							std::cout << "found " << attribute << fishnames[z - 2] << " of hint " << (h + 1) << " in " << sets;
							handle_fish(s, results[r].positions, (1ul << h), changed);
							std::cout << std::endl;
							if(changed)
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
			positions_t          urrows;
			positions_t          urclms;
			candidates_t         urcandidates;
			int                  urtype;
			std::vector<index_t> ercells;
			candidates_t         ercandidates;
		};

		// Finds all cells that have more candidates than the other three.
		positions_t excessCells(candidates_t h[4], candidates_t& leastCommon) const
		{
			positions_t result = 0;
			leastCommon = h[0] & h[1] & h[2] & h[3];
			for (int i = 0; i < 4; i++)
			{
				if (h[i] & ~leastCommon)
					result |= (1ul << i);
			}
			return result;
		}
		void find_ur(Sudoku* s,std::vector<urResult>& results) const
		{
			auto p0=positions_t(1);

			index_t ur[4];
			for(index_t i = 0; i < s->nn(); i++)
			{
				auto indexerR1 = s->getIndexer(Sudoku::Indexer::Type::row,i);
				for(index_t ii = i+1; ii < s->nn(); ii++)
				{
					auto indexerR2 = s->getIndexer(Sudoku::Indexer::Type::row,ii);
					for(index_t j = 0; j < s->nn(); j++)
					{
						ur[0] = (*indexerR1)(j);
						ur[1] = (*indexerR2)(j);
						for(index_t jj = j + 1; jj < s->nn(); jj++)
						{
							ur[2] = (*indexerR1)(jj);
							ur[3] = (*indexerR2)(jj);
							if(bitcount((1ul << s->getBox(ur[0])) |
										(1ul << s->getBox(ur[1])) |
										(1ul << s->getBox(ur[2])) |
										(1ul << s->getBox(ur[3]))) == 2) // only two boxed covered?
							{
								candidates_t candidates[4], urcandidates;
								for(int k = 0; k < 4; k++)
									candidates[k] = s->cells[ur[k]].candidates;

								positions_t extraPositions = excessCells(candidates, urcandidates);

								if(bitcount(urcandidates) == 2){
									switch(extraPositions)
									{
									case 1:
									case 2:
									case 4:
									case 8:
										// type 1
										results.push_back({static_cast<positions_t>((p0 << i) | (p0 << ii)),
														   static_cast<positions_t>((p0 << j) | (p0 << jj)),
														   urcandidates, 1, {ur[get_last_set(extraPositions)]}, urcandidates});
										break;
									case 3:
									case 5:
									case 10:
									case 12:
										// type 2
										auto exx = bitindex<4>(candidates, extraPositions);
										if(bitcount(bitunion(exx) & ~urcandidates) == 1) // only one extra candidate?
										{
											if(candidates_t extraCandidates = bitintersection(exx) & ~urcandidates)
											{
												std::vector<index_t> delpos;
												visible(s, bitindex<4>(ur, extraPositions), delpos, false);
												results.push_back({static_cast<positions_t>((1ul << i) | (1ul << ii)),
																   static_cast<positions_t>((1ul << j) | (1ul << jj)),
																   urcandidates, 2, delpos, extraCandidates});
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
		}
		void handle_ur(Sudoku* s, urResult urr, bool& changed) const
		{
			// highlight the board
			for (int i = 0; i < s->nn(); i++)
				if (urr.urrows & (1ul << i))
				{
					s->markers[{-1, i, 0}] = { 1 };
					for (int j = 0; j < s->nn(); j++)
						if (urr.urclms & (1ul << j))
						{
							s->markers[{j, -1, 0}] = { 1 };
							for (int h = 0; h < s->nn(); h++)
								if (urr.urcandidates & (1ul << h))
								{
									s->markers[{j, i, h + 1}] = { 0 };
								}
						}
				}

			// erase
			int nelims = 0;
			for (auto c : urr.ercells)
			{
				uint16_t& cellcandidates = s->cells[c].candidates;
				if (!cellcandidates)
					continue;
				nelims += bitcount(cellcandidates);
				cellcandidates &= ~urr.ercandidates;
				nelims -= bitcount(cellcandidates);
				for (int h = 0; h < s->nn(); h++)
					if (urr.ercandidates & (1ul << h))
						s->markers[{s->getClm(c), s->getRow(c), h + 1}] = { 2 };
			}
			std::cout << " -> eliminate " << nelims << " occurences";
			for (auto c : urr.ercells)
			{
				// set
				uint16_t& cellcandidates = s->cells[c].candidates; 
				if (!cellcandidates)
					continue;
				if ((cellcandidates & (cellcandidates - 1)) == 0) /// only one bit set?
				{
					std::cout << " -> set " << bitsToString(cellcandidates, 1) << " at r" << s->getRow(c) + 1ul << "c" << s->getClm(c) + 1;
					s->markers[{s->getClm(c), s->getRow(c), 0}] = { 0 };

					s->cells[c].value = cellcandidates;
					eliminateCandidates(s, cellcandidates,
										{s->getIndexer(Sudoku::Indexer::Type::row, s->getRow(c)),
										 s->getIndexer(Sudoku::Indexer::Type::clm, s->getClm(c)),
										 s->getIndexer(Sudoku::Indexer::Type::box, s->getBox(c))});
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
				std::cout << "found Unique Rectangle Type " << a.urtype << " of candidates " << bitsToString(a.urcandidates, 1) << " in " 
					<< descr(Sudoku::Indexer::Type::row) << bitsToString(a.urrows, 1, "")
					<< descr(Sudoku::Indexer::Type::clm) << bitsToString(a.urclms, 1, "");

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
		bool brut(Sudoku* s, Sudoku s0, int startcell) const
		{
			//iters++;
			//if (startcell == s->nn() * s->nn())
			//{
			//	// all cells successfully filled
			//	*s=s0;
			//	return true;
			//}
			//if (s0.cells[startcell].value) return brut(s, s0, startcell + 1);
			//for (index_t h = 0; h < s->nn(); h++)
			//{
			//	if (s0.cells[startcell].candidates & (1ul << h))
			//	{
			//		Sudoku s1 = s0;
			//		s1.cells[startcell].value = (1ul << h);
			//		eliminateCandidates(&s1, (1ul << h),
			//							{s1.getIndexer(Sudoku::Indexer::Type::row, s1.getRow(startcell)),
			//							 s1.getIndexer(Sudoku::Indexer::Type::clm, s1.getClm(startcell)),
			//							 s1.getIndexer(Sudoku::Indexer::Type::box, s1.getBox(startcell))});

			//		if (brut(s, s1, startcell + 1))
			//			return true;
			//	}
			//}
			return false;
		}

	public:
		void apply(Sudoku* s, int maxorder) const override
		{
			if (brut(s, *s, 0))
				std::cout << setColor(lgreen) << "Solution found! Iterations: " << iters << setColor() << std::endl;
			else
				std::cout << setColor(lred) << "There is no solution! Iterations: " << iters << setColor() << std::endl;
		}
	};
}
