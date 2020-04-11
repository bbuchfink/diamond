#pragma once
#include <unordered_map>
#include <vector>

template<typename T>
class SimpleMatrix {
public:
	virtual T get_max_elm() = 0;
	virtual inline T get_elm(uint32_t i_idx, uint32_t j_idx) = 0;
	virtual inline void set_elm(uint32_t i_idx, uint32_t j_idx, T element) = 0;
	virtual T& at(uint32_t i_idx, uint32_t j_idx) = 0;
	virtual T& at(uint64_t c_idx) = 0;
	virtual inline uint64_t get_nrows() = 0;
	virtual inline uint64_t get_ncols() = 0;
	virtual SimpleMatrix<T>* get_new_instance(int nrows, int ncols) = 0;
	virtual SimpleMatrix<T>* get_new_instance() = 0;
	virtual void print(){
		for( uint32_t irow=0; irow < get_nrows(); irow++){
			for( uint32_t icol=0; icol < get_ncols(); icol++){
				printf(" %.6f",get_elm(irow,icol));
			}
			printf("\n");
		}
	}

	virtual inline uint64_t get_idx(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < get_nrows());
		assert(j_idx < get_ncols());
		return i_idx * get_ncols() + j_idx;
	}

	virtual inline void get_indices (uint64_t combined, uint32_t* i, uint32_t* j){
		*j = combined % get_ncols();
		*i = ( combined - *j) / get_ncols();
	}

	virtual SimpleMatrix<T>* minus(SimpleMatrix<T>* m) = 0;
	virtual SimpleMatrix<T>* base_minus(SimpleMatrix<T>* m){
		assert(this->get_nrows() == m->get_nrows());
		assert(this->get_ncols() == m->get_ncols());
		SimpleMatrix<T>* result = this->get_new_instance(this->get_nrows(),m->get_ncols());
		for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
			for(uint32_t icol = 0; icol<m->get_ncols(); icol++){
				result->set_elm(irow, icol, this->get_elm(irow, icol) - m->get_elm(irow,icol));
			}
		}
		return result;
	}
	virtual SimpleMatrix<T>* plus(SimpleMatrix<T>* m) = 0;
	virtual SimpleMatrix<T>* base_plus(SimpleMatrix<T>* m){
		assert(this->get_nrows() == m->get_nrows());
		assert(this->get_ncols() == m->get_ncols());
		SimpleMatrix<T>* result = this->get_new_instance(this->get_nrows(),m->get_ncols());
		for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
			for(uint32_t icol = 0; icol<m->get_ncols(); icol++){
				result->set_elm(irow, icol, this->get_elm(irow, icol) + m->get_elm(irow, icol));
			}
		}
		return result;
	}
	virtual SimpleMatrix<T>* multiply(SimpleMatrix<T>* m) = 0;
	virtual SimpleMatrix<T>* base_multiply(SimpleMatrix<T>* m){
		printf("standard mult\n");
		SimpleMatrix<T>* result = this->get_new_instance(this->get_nrows(),m->get_ncols());
		for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
			for(uint32_t icol = 0; icol<m->get_ncols(); icol++){
				result->set_elm(irow, icol, 0);
			}
		}
		for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
			for(uint32_t icol = 0; icol<m->get_ncols(); icol++){
				for(uint32_t k = 0; k<m->get_nrows(); k++){
					result->set_elm(irow, icol, result->get_elm(irow, icol) + this->get_elm(irow, k) * m->get_elm(k,icol));
				}
			}
		}
		return result;
	}
	virtual double norm(double p, double q){
		double result = 0.0;
		double exp1 = q/p;
		for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
			double colnorm = 0.0;
			for(uint32_t icol = 0; icol<this->get_ncols(); icol++){
					colnorm += std::pow(std::abs(this->get_elm(irow, icol)), p);
			}
			result += std::pow(colnorm, exp1);
		}
		return std::pow(result, 1.0/q);
	}

	virtual ~SimpleMatrix() {};
};

template<typename T>
class SymmetricSimpleMatrix : virtual public SimpleMatrix<T> {
public:
	static inline uint64_t get_lower_index (uint32_t i_idx, uint32_t j_idx){
		uint32_t min = std::min(i_idx,j_idx);
		uint32_t max = std::max(i_idx,j_idx);
		return  max*((uint64_t) max+1)/2 + min;
	}
	static inline void get_symmetric_index (uint64_t combined, uint32_t* max, uint32_t* min){
		// solving max^2-max-2*combined = 0. -1 to convert to indexing varaible in C, avoiding overflow
		// by removing an additional factor of 2 from under the square root
		*max = (uint32_t) std::floor(0.5+2*sqrt(0.0625+0.5*combined)) - 1;
		*min = combined - (*max*((uint64_t)*max+1))/2;
	}
};

template<typename T>
class SparseSimpleMatrix : virtual public SimpleMatrix<T> {
public:
	virtual uint64_t get_n_nonzero_elements() = 0;
	virtual inline bool has_elm(uint32_t i_idx, uint32_t j_idx) = 0;
	virtual typename std::unordered_map<uint64_t, T>::iterator begin() = 0;
	virtual typename std::unordered_map<uint64_t, T>::iterator end() = 0;
};

template<typename T>
class SparseSimpleMatrixImpl : public SparseSimpleMatrix<T> {
private:
	uint32_t nrows, ncols;
	T tol;
	std::unordered_map<uint64_t,T> mat;

public:
	SparseSimpleMatrixImpl(int nrows, int ncols) : nrows(nrows), ncols(ncols), tol((T) 1e-12) {};
	SparseSimpleMatrixImpl(int nrows, int ncols, T tol) : nrows(nrows), ncols(ncols), tol(std::abs(tol)) {};

	virtual SimpleMatrix<T>*  get_new_instance(int nrows, int ncols){
		return new SparseSimpleMatrixImpl<T>(nrows, ncols);
	}
	virtual SimpleMatrix<T>*  get_new_instance(){
		return new SparseSimpleMatrixImpl<T>(nrows, ncols);
	}

	virtual ~SparseSimpleMatrixImpl(){ mat.clear(); };

	virtual inline void set_elm(uint32_t i_idx, uint32_t j_idx, T element){
		assert(i_idx < nrows);
		assert(j_idx < ncols);
		if(std::abs(element) > tol) {
			mat[SimpleMatrix<T>::get_idx(i_idx,j_idx)] = element;
		}
	}
	virtual inline T get_elm(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < nrows);
		assert(j_idx < ncols);
		auto it = mat.find(SimpleMatrix<T>::get_idx(i_idx,j_idx));
		return it != mat.end() ? it-> second : (T) 0.0;
	}
	virtual inline bool has_elm(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < nrows);
		assert(j_idx < ncols);
		auto it = mat.find(SimpleMatrix<T>::get_idx(i_idx,j_idx));
		return it != mat.end();
	}
	virtual inline uint64_t get_nrows(){
		return nrows;
	};
	virtual inline uint64_t get_ncols(){
		return ncols;
	};
	virtual inline uint64_t get_n_nonzero_elements(){
		return mat.size();
	};
	virtual typename std::unordered_map<uint64_t,T>::iterator begin(){
	 	return mat.begin();
	};
	virtual typename std::unordered_map<uint64_t,T>::iterator end(){
		return mat.end();
	};
	virtual T& at(uint32_t i_idx, uint32_t j_idx){
		return at(SimpleMatrix<T>::get_idx(i_idx,j_idx));
	};
	virtual T& at(uint64_t c_idx){
		assert(c_idx < nrows*ncols);
		return mat[c_idx];
	};

	void purge(){
		auto it = mat.begin();
		while(it!=mat.end()){
			if(std::abs(it->second) < tol){
				it = mat.erase(it);
			}
			else{
				it++;
			}
		}
	}
	virtual T get_max_elm(){
		T m = 0.0;
		for(auto& it : mat){
			m = std::max(m, it.second);
		}
		return m;
	};
	virtual SimpleMatrix<T>* minus(SimpleMatrix<T>* m){
		if(SparseSimpleMatrix<T>* n = dynamic_cast<SparseSimpleMatrix<T>*> (m)){
			// TODO: Imporove this implementation
			SparseSimpleMatrix<T>* result = new SparseSimpleMatrixImpl(this->get_nrows(),n->get_ncols());

			std::vector<std::set<uint32_t>> k_to_j(n->get_nrows());
			for(auto it_m = n->begin(); it_m != n->end(); it_m++){
				uint32_t irow,icol;
				n->get_indices(it_m->first, &irow, &icol);
				result->set_elm(irow, icol, this->get_elm(irow,icol) - it_m->second);
			}
			for(auto it_this = this->begin(); it_this != this->end(); it_this++){
				uint32_t irow,icol;
				this->get_indices(it_this->first, &irow, &icol);
				if(!n->has_elm(irow,icol)) result->set_elm(irow, icol, it_this->second);
			}
			return result;
		}
		return this->base_minus(m);
	}
	virtual SimpleMatrix<T>* plus(SimpleMatrix<T>* m){
		if(SparseSimpleMatrix<T>* n = dynamic_cast<SparseSimpleMatrix<T>*> (m)){
			// TODO: Imporove this implementation
			SparseSimpleMatrix<T>* result = new SparseSimpleMatrixImpl(this->get_nrows(),n->get_ncols());

			std::vector<std::set<uint32_t>> k_to_j(n->get_nrows());
			for(auto it_m = n->begin(); it_m != n->end(); it_m++){
				uint32_t irow,icol;
				n->get_indices(it_m->first, &irow, &icol);
				result->set_elm(irow, icol, this->get_elm(irow,icol) + it_m->second);
			}
			for(auto it_this = this->begin(); it_this != this->end(); it_this++){
				uint32_t irow,icol;
				this->get_indices(it_this->first, &irow, &icol);
				if(!n->has_elm(irow,icol)) result->set_elm(irow, icol, it_this->second);
			}
			return result;
		}
		return this->base_plus(m);
	}
	virtual SimpleMatrix<T>* multiply(SimpleMatrix<T>* m){
		if(SparseSimpleMatrix<T>* n = dynamic_cast<SparseSimpleMatrix<T>*> (m)){
			SparseSimpleMatrix<T>* result = new SparseSimpleMatrixImpl(this->get_nrows(),n->get_ncols());
			std::vector<std::vector<uint32_t>> k_to_j(n->get_nrows());
			for(uint32_t k = 0; k<n->get_nrows(); k++){
				k_to_j[k] = std::vector<uint32_t>();
			}
			for(auto it_m = n->begin(); it_m != n->end(); it_m++){
				uint32_t k2,icol;
				n->get_indices(it_m->first, &k2, &icol);
				k_to_j[k2].push_back(icol);
			}

			auto it_this = this->begin();
			uint32_t ithis = 0;
			while(it_this != this->end()){
				uint32_t irow, k;
				this->get_indices(it_this->first, &irow, &k);
				for( uint32_t const & icol : k_to_j[k] ){
					T el = result->get_elm(irow, icol);
					result->set_elm(irow, icol, el + it_this->second * n->get_elm(k, icol));
				}
				it_this++;
			}
			return result;
		}
		return this->base_multiply(m);
	}
};

template<typename T>
class DenseSimpleMatrix : public SimpleMatrix<T> {
private:
	uint32_t nrows, ncols;
	T* mat;
public:
	DenseSimpleMatrix(int nrows, int ncols) : nrows(nrows), ncols(ncols) {
		mat = new T[nrows * ncols];
	};
	DenseSimpleMatrix(int nrows, int ncols, T init) : nrows(nrows), ncols(ncols) {
		mat = new T[nrows * ncols];
		for(uint64_t i = 0; i<nrows*ncols; i++){
			mat[i] = init;
		}
	};
	virtual ~DenseSimpleMatrix(){
		delete [] mat;
		mat = nullptr;
	};
	virtual SimpleMatrix<T>* get_new_instance(int nrows, int ncols){
		return new DenseSimpleMatrix<T>(nrows, ncols);
	}
	virtual SimpleMatrix<T>* get_new_instance(){
		return new DenseSimpleMatrix<T>(nrows, ncols);
	}

	virtual inline void set_elm(uint32_t i_idx, uint32_t j_idx, T element){
		mat[SimpleMatrix<T>::get_idx(i_idx,j_idx)] = element;
	}
	virtual inline T get_elm(uint32_t i_idx, uint32_t j_idx){
		return mat[SimpleMatrix<T>::get_idx(i_idx,j_idx)];
	}
	virtual inline uint64_t get_nrows(){
		return nrows;
	};
	virtual inline uint64_t get_ncols(){
		return ncols;
	};
	virtual T& at(uint32_t i_idx, uint32_t j_idx){
		return at(SimpleMatrix<T>::get_idx(i_idx,j_idx));
	};
	virtual T& at(uint64_t c_idx){
		assert(c_idx < nrows * ncols);
		return mat[c_idx];
	};
	virtual T get_max_elm(){
		// TODO: implement me
		return 0;//*std::max_element(mat,mat+n_max);
	};

	SimpleMatrix<T>* multiply(SimpleMatrix<T>* m){
		if(SparseSimpleMatrix<T>* n = dynamic_cast<SparseSimpleMatrix<T>*> (m)){
			DenseSimpleMatrix<T>* result = new DenseSimpleMatrix(this->get_nrows(),n->get_ncols());
			for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
				for(uint32_t icol = 0; icol<n->get_ncols(); icol++){
					result->set_elm(irow, icol, 0);
				}
			}
			for(uint32_t irow = 0; irow<this->get_nrows(); irow++){
				auto it = n->begin();
				while(it!=n->end()){
					uint32_t k,icol;
					n->get_indices (it->first, &k, &icol);
					result->at(irow, icol) += this->get_elm(irow, k) * it->second;
					it++;
				}
			}
			return result;
		}
		return this->base_multiply(m);
	}
};

template<typename T>
class SparseSymmetricSimpleMatrix : public SparseSimpleMatrix<T>, public SymmetricSimpleMatrix<T> {
private:
	uint32_t n_max;
	T tol;
	std::unordered_map<uint64_t,T> mat;

public:
	SparseSymmetricSimpleMatrix(int n) : n_max(n), tol((T) 1e-12) {};
	SparseSymmetricSimpleMatrix(int n, T tol) : n_max(n), tol(std::abs(tol)) {};

	virtual ~SparseSymmetricSimpleMatrix(){ mat.clear(); };
	virtual SimpleMatrix<T>* get_new_instance(int nrows, int ncols){
		return new SparseSymmetricSimpleMatrix<T>(nrows, tol);
	}
	virtual SimpleMatrix<T>* get_new_instance(){
		return new SparseSymmetricSimpleMatrix<T>(n_max, tol);
	}

	virtual inline void set_elm(uint32_t i_idx, uint32_t j_idx, T element){
		if( std::abs(element) > tol ) {
			assert(i_idx < n_max);
			assert(j_idx < n_max);
			const uint64_t my_idx=this->get_lower_index(i_idx,j_idx);
			mat[my_idx] = element;
		}
	}
	virtual inline T get_elm(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t my_idx=this->get_lower_index(i_idx,j_idx);
		auto it = mat.find(my_idx);
		return it != mat.end() ? it-> second : (T) 0.0;
	}
	virtual inline bool has_elm(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t my_idx=this->get_lower_index(i_idx,j_idx);
		auto it = mat.find(my_idx);
		return it != mat.end();
	}
	virtual inline uint64_t get_nrows(){
		// TODO: this can be determined from mat (since symmetric) then the matrix can be unbounded!
		return n_max;
	};
	virtual inline uint64_t get_ncols(){
		// TODO: this can be determined from mat (since symmetric) then the matrix can be unbounded!
		return n_max;
	};
	virtual inline uint64_t get_n_nonzero_elements(){
		return mat.size();
	};
	virtual typename std::unordered_map<uint64_t,T>::iterator begin(){
		return mat.begin();
	};
	virtual typename std::unordered_map<uint64_t,T>::iterator end(){
		return mat.end();
	};
	virtual T& at(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t c_idx=this->get_lower_index(i_idx,j_idx);
		return at(c_idx);
	};
	virtual T& at(uint64_t c_idx){
		assert(c_idx < n_max*(n_max+1)/2);
		return mat[c_idx];
	};

	void purge(){
		auto it = mat.begin();
		while(it!=mat.end()){
			if (std::abs(it->second) < tol ){
				it = mat.erase(it);
			}
			else{
				it++;
			}
		}
	}

	virtual T get_max_elm(){
		T m = 0.0;
		for ( auto& it : mat ) {
			m = std::max(m,it.second);
		}
		return m;
	};

	void resize(uint32_t new_n_max){
		n_max = new_n_max;
	}
	// TODO: implement valid matrix ops
	virtual SimpleMatrix<T>* minus(SimpleMatrix<T>* m){
		return this->base_minus(m);
	}
	virtual SimpleMatrix<T>* plus(SimpleMatrix<T>* m){
		return this->base_plus(m);
	}
	virtual SimpleMatrix<T>* multiply(SimpleMatrix<T>* m){
		return this->base_multiply(m);
	}
};

template<typename T>
class DenseSymmetricSimpleMatrix : public SymmetricSimpleMatrix<T> {
private:
	uint32_t n_max;
	T* mat;
public:
	DenseSymmetricSimpleMatrix(int n) : n_max(n) {
		mat = new T[(n_max*(n_max+1))/2];
	};
	virtual ~DenseSymmetricSimpleMatrix(){
		delete [] mat;
		mat = nullptr;
	};
	virtual SimpleMatrix<T>* get_new_instance(int nrows, int ncols){
		return new DenseSymmetricSimpleMatrix<T>(nrows * ncols);
	}
	virtual SimpleMatrix<T>* get_new_instance(){
		return new DenseSymmetricSimpleMatrix<T>(n_max);
	}

	virtual inline void set_elm(uint32_t i_idx, uint32_t j_idx, T element){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t my_idx=SymmetricSimpleMatrix<T>::get_lower_index(i_idx,j_idx);
		mat[my_idx] = element;
	}
	virtual inline T get_elm(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t my_idx=SymmetricSimpleMatrix<T>::get_lower_index(i_idx,j_idx);
		return mat[my_idx];
	}
	virtual inline uint64_t get_nrows(){
		return n_max;
	};
	virtual inline uint64_t get_ncols(){
		return n_max;
	};
	virtual T& at(uint32_t i_idx, uint32_t j_idx){
		assert(i_idx < n_max);
		assert(j_idx < n_max);
		const uint64_t c_idx=SymmetricSimpleMatrix<T>::get_lower_index(i_idx,j_idx);
		return at(c_idx);
	};
	virtual T& at(uint64_t c_idx){
		assert(c_idx < n_max*(n_max+1)/2);
		return mat[c_idx];
	};
	virtual T get_max_elm(){
		return *std::max_element(mat,mat+n_max);
	};
};
