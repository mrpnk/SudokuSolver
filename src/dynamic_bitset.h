

template<int block>
class dynamic_bitset{

	std::vector<bool> v;
public:

	dynamic_bitset(int n = 0,bool value = false){
		resize(n,value);
	}

	void resize(size_t n, bool value = false){
		v.resize(n,value);
	}

	void fill(bool value){
		v=std::vector<bool>(v.size(), value);
	}

	void set(size_t i, bool value = true){
		v[i]=value;
	}

	void and_assign(size_t i, bool value){
		v[i]=v[i]&&value;
	}


	bool operator[](size_t i) const { return v[i]; }

	size_t size() const {
		return v.size();
	}

	bool any() const {
		for (bool b : v)
			if (b)
				return true;
		return false;
	}

	int popcount() const {
		int count{ 0 };
		for (const auto& a : v)
			if (a) count++;

		return count;
	}

	std::string to_string() const{
		std::string s(v.size(),'.');
		for(int i=0;bool b : v){
			if(b)s[i]='X';
			++i;
		}
		return s;
	}
};


using dyn_bs = dynamic_bitset<1>;
