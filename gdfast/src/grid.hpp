#pragma once


class Grid1dIrregular {
	public:
		Grid1dIrregular(double_vector rborders) : rborders(rborders), index(0) {
			rbordersp = rborders.data().begin();
			_length = rborders.size();
		}
		int findindex(double r) {
			if( (r >= rbordersp[index]) && (r < rbordersp[index+1]))
				return index;
			if (r < rbordersp[0])
				return 0;
			if (r >= rbordersp[_length-1])
				return _length-1;
	
			if(r < rbordersp[index]) // we should go to the left 
				while(r < rbordersp[index]) // if r >= rbordersp[index], we found it
					index--;
			else
				while(r >= rbordersp[index+1]) // if r < rbordersp[index+1], we found it
					index++;
			return index;
		}
		bool inrange(double r) {
			return (r >= rbordersp[0]) && (r < rbordersp[_length-1]); 
		}
		int length() { return _length; }
	private:
		int index; // remember the old position (this is slow for random r, but slow for slowly changing r), index refers to the 'left' edge
		int _length;
		double_vector rborders;
		double* rbordersp;
};

class Grid1dRegular {
	public:
		Grid1dRegular(double x1, double x2, int length) : x1(x1), x2(x2), _length(length) {}
		int findindex(double x) {
			
			return x < x1 ? 0 :
					(x >= x2 ? _length-1 : (int)((x-x1)/(x2-x1)*(_length)));
		}
		bool inrange(double r) {
			return (r >= x1) && (r < x2); 
		}
		int border(int i) {
			return x1 + (x2-x1)*i/length; 
		}
		int length() { return length; }
		double x1, x2;
		int length;
};

template<int D, int N>
class GridLinear {
	
};



/*
class Grid1dRegular  {
	public:
		Grid1dRegular(double x1, double x2, int length) : x1(x1), x2(x2), _length(length) {}
		virtual int findindex(double x) {
			
			return x < x1 ? 0 :
					(x >= x2 ? _length-1 : (int)((x-x1)/(x2-x1)*(_length)));
		}
		virtual bool inrange(double r) {
			return (r >= x1) && (r < x2); 
		}
		virtual int border(int i) {
			return x1 + (x2-x1)*i/length; 
		}
		virtual int length() { return length; }
		double x1, x2;
		int length;
};

class Grid1dRegularLinear : public Grid1dRegular {
	public:
		Grid1dRegular(double x1, double x2, int length) : Grid1dRegular(x1, x2, length)
};

class Grid1dRegularLog : public Grid1dRegular {
	public:
		Grid1dRegular(double x1, double x2, int length) : Grid1dRegular(x1, x2, length) {}
		virtual double f(double x) {
			return pow(10, x);
		}
		virtual double f_inv(double u) {
			return log10(u);
		}
		virtual double df_invdx(doube u) {
			return 1/(log(10)*u);
		}
};
template<class Base=Aperture1d>
class Grid1dRegularLog : public Base {
	public:
		Grid1dRegularLog(double u1, double u2, int length) : u1(u1), u2(u2), _length(length) {}
		int findindex(double x) {
			double u = log10(x);
			return u < u1 ? 0 :
					(u >= u2 ? -1 : (int)((u-u1)/(u2-u1)*(_length)));
		}
		bool inrange(double x) {
			double u = log10(x);
			return (u >= u1) && (u < u2); 
		}
		double uniform_transform(double uniform) {
			return pow(10, uniform * (u2-u1) + u1);
		}
		int length() { return _length; }
		double u1, u2;
		int _length;
};

template<int D, int N>
class GridLinear {
	
}; * */