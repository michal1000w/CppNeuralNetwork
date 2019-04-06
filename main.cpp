#include <iostream>
#include <cmath>
#include <conio.h>
#include <ctime>
#include <vector>
#include <string>


using namespace std;

class Matrix {
public:
	friend class NeuralNetwork;
	friend class NeuNet;
	friend Matrix* operator += (const Matrix* lhs, const Matrix& rhs);
	Matrix(unsigned int, unsigned int, double**);
	Matrix(unsigned int, unsigned int);
	Matrix(unsigned int);
	Matrix();
	void add(const double*, ...);
	void add(unsigned int, ...);
	void add(string);
	Matrix print();
	Matrix print(short);
	Matrix T();
	Matrix* t();

	double** getArray();

	//funkcjonalnoœæ pod AI
	Matrix expa(bool inverted = 0);
	Matrix sigmoid();
	Matrix sigmoid_derivative();
	//

	Matrix operator + (const Matrix& rhs);
	Matrix operator + (const Matrix* rhs);
	Matrix operator += (const Matrix* rhs);
	Matrix operator += (const Matrix& rhs);
	Matrix operator - (const Matrix* rhs);
	Matrix operator - (const Matrix& rhs);
	Matrix operator * (const Matrix& rhs);
	Matrix operator * (const Matrix* rhs);
	Matrix operator * (const double& rhs);
	Matrix operator *= (const Matrix& rhs);
	void operator delete(void* ptr);
private:
	void initMatrix();
	unsigned int liczba_macierzy;
	unsigned int liczba_elementow;
	double** arrays;
};

void Matrix::operator delete(void* ptr) {
	delete (ptr);
}

double** Matrix::getArray() {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++) newArray[i] = new double[x];

	for (int j = 0; j < y; j++) {
		for (int i = 0; i < x; i++) {
			newArray[j][i] = this->arrays[j][i];
		}
	}

	return newArray;
}

Matrix Matrix::expa(bool inverted) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++) newArray[i] = new double[x];

	if (inverted)
		for (int j = 0; j < y; j++)
			for (int i = 0; i < x; i++)
				newArray[j][i] = exp((this->arrays[j][i]) * -1);
	else
		for (int j = 0; j < y; j++)
			for (int i = 0; i < x; i++)
				newArray[j][i] = exp(this->arrays[j][i]);

	return Matrix(y, x, newArray);
}

Matrix Matrix::sigmoid() {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++) newArray[i] = new double[x];

	newArray = this->expa(1).arrays;

	for (int j = 0; j < y; j++) {
		for (int i = 0; i < x; i++) {
			newArray[j][i] = 1 / (1 + newArray[j][i]);
		}
	}

	return Matrix(y, x, newArray);
}

Matrix Matrix::sigmoid_derivative() {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++) newArray[i] = new double[x];

	for (int j = 0; j < y; j++) {
		for (int i = 0; i < x; i++) {
			newArray[j][i] = this->arrays[j][i] * (1 - this->arrays[j][i]);
		}
	}

	return Matrix(y, x, newArray);
}

Matrix Matrix::T() {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[x];
	for (int i = 0; i < x; i++) {
		newArray[i] = new double[y];
	}

	for (int j = 0; j < y; j++) {
		for (int i = 0; i < x; i++) {
			newArray[i][j] = this->arrays[j][i];
		}
	}

	return Matrix(x, y, newArray);
}

Matrix* Matrix::t() {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[x];
	for (int i = 0; i < x; i++) {
		newArray[i] = new double[y];
	}

	for (int j = 0; j < y; j++) {
		for (int i = 0; i < x; i++) {
			newArray[i][j] = this->arrays[j][i];
		}
	}

	Matrix* output = new Matrix(x, y, newArray);


	return output;
}

Matrix* operator += (const Matrix* lhs, const Matrix& rhs) {
	if (!(lhs->liczba_elementow == rhs.liczba_elementow && lhs->liczba_macierzy == rhs.liczba_macierzy)) cout << "Nie mozna dodac tych macierzy" << endl;
	else {

		int y = lhs->liczba_macierzy;
		int x = lhs->liczba_elementow;

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				lhs->arrays[i][j] = lhs->arrays[i][j] + rhs.arrays[i][j];
			}
		}
	}

	Matrix* output = new Matrix(lhs->liczba_macierzy, lhs->liczba_elementow, lhs->arrays);

	return output;
}

Matrix Matrix::operator +(const Matrix& rhs) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++)
		newArray[i] = new double[x];

	if (!(this->liczba_elementow == rhs.liczba_elementow && this->liczba_macierzy == rhs.liczba_macierzy)) cout << "Nie mozna dodac tych macierzy" << endl;
	else {

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				newArray[i][j] = this->arrays[i][j] + rhs.arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, newArray);
}

Matrix Matrix::operator +(const Matrix* rhs) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++)
		newArray[i] = new double[x];

	if (!(this->liczba_elementow == rhs->liczba_elementow && this->liczba_macierzy == rhs->liczba_macierzy)) cout << "Nie mozna dodac tych macierzy" << endl;
	else {

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				newArray[i][j] = this->arrays[i][j] + rhs->arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, newArray);
}

Matrix Matrix::operator +=(const Matrix* rhs) {
	if (!(this->liczba_elementow == rhs->liczba_elementow && this->liczba_macierzy == rhs->liczba_macierzy)) cout << "Nie mozna dodac tych macierzy" << endl;
	else {

		int y = this->liczba_macierzy;
		int x = this->liczba_elementow;

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				this->arrays[i][j] = this->arrays[i][j] + rhs->arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
}

Matrix Matrix::operator +=(const Matrix& rhs) {
	if (!(this->liczba_elementow == rhs.liczba_elementow && this->liczba_macierzy == rhs.liczba_macierzy)) cout << "Nie mozna dodac tych macierzy" << endl;
	else {

		int y = this->liczba_macierzy;
		int x = this->liczba_elementow;

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				this->arrays[i][j] = this->arrays[i][j] + rhs.arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
}


Matrix Matrix::operator -(const Matrix& rhs) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++)
		newArray[i] = new double[x];

	if (!(this->liczba_elementow == rhs.liczba_elementow && this->liczba_macierzy == rhs.liczba_macierzy)) cout << "Nie mozna odjac tych macierzy" << endl;
	else {

		for (unsigned int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				newArray[i][j] = this->arrays[i][j] - rhs.arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, newArray);
}

Matrix Matrix::operator -(const Matrix* rhs) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++)
		newArray[i] = new double[x];

	if (!(this->liczba_elementow == rhs->liczba_elementow && this->liczba_macierzy == rhs->liczba_macierzy)) cout << "Nie mozna odjac tych macierzy" << endl;
	else {

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				newArray[i][j] = this->arrays[i][j] - rhs->arrays[i][j];
			}
		}
	}

	return Matrix(this->liczba_macierzy, this->liczba_elementow, newArray);
}

Matrix Matrix::operator *=(const Matrix& rhs) {
	unsigned int y = this->liczba_macierzy;
	unsigned int x = this->liczba_elementow;

	double** newArray = new double*[y];
	for (int i = 0; i < y; i++)
		newArray[i] = new double[x];

	if (!(this->liczba_elementow == rhs.liczba_elementow && this->liczba_macierzy == rhs.liczba_macierzy)) cout << "Nie mozna zrobic *=" << endl;
	else {

		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				newArray[i][j] = this->arrays[i][j] * rhs.arrays[i][j];
			}
		}
	}

	//czyszczenie
	for (int i = 0; i < y; i++) delete[] rhs.arrays[i];
	for (int i = 0; i < this->liczba_macierzy; i++) delete[] this->arrays[i];
	delete[] this->arrays;
	delete[] rhs.arrays;

	return Matrix(this->liczba_macierzy, this->liczba_elementow, newArray);
}

Matrix Matrix::operator *(const Matrix& rhs) {
	if (!(this->liczba_elementow == rhs.liczba_macierzy)) {
		cout << "Nie mozna pomnozyc tych macierzy" << endl;
		return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
	}
	else {
		int y1 = this->liczba_macierzy;
		int x1 = this->liczba_elementow;
		int y2 = rhs.liczba_macierzy;
		int x2 = rhs.liczba_elementow;

		double** newArray = new double*[y1];
		for (int i = 0; i < y1; i++) {
			newArray[i] = new double[x2];
		}

		double suma = 0;
		int w, k, t;
		w = k = 0;
		t = rhs.liczba_elementow * this->liczba_macierzy;
		int t1 = t;

		while (t > 0) {
			suma = 0;
			for (int i = 0; i < x1; i++) {
				suma += this->arrays[w][i] * rhs.arrays[i][k];
			}

			newArray[w][k] = suma;

			k = (k + 1)%x2;
			t--;
			if (t%x2 == 0 && t != t1 && w < (y1-1)) w++;
			//t--;
			if (t == 0) break;
		}

		//czyszczenie
		/*for (int i = 0; i < y2; i++) delete[] rhs.arrays[i];
		for (int i = 0; i < y1; i++) delete[] this->arrays[i];
		delete[] this->arrays;
		delete[] rhs.arrays;*/

		return Matrix(y1, x2, newArray);
	}
}

Matrix Matrix::operator *(const Matrix* rhs) {
	if (!(this->liczba_elementow == rhs->liczba_macierzy)) {
		cout << "Nie mozna pomnozyc tych macierzy" << endl;
		return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
	}
	else {
		int y1 = this->liczba_macierzy;
		int x1 = this->liczba_elementow;
		int y2 = rhs->liczba_macierzy;
		int x2 = rhs->liczba_elementow;

		double** newArray = new double*[y1];
		for (int i = 0; i < y1; i++) {
			newArray[i] = new double[x2];
		}

		double suma = 0;
		int w, k, t;
		w = k = 0;
		t = rhs->liczba_elementow * this->liczba_macierzy;
		int t1 = t;

		while (t > 0) {
			suma = 0;
			for (int i = 0; i < x1; i++) {
				suma += this->arrays[w][i] * rhs->arrays[i][k];
			}

			newArray[w][k] = suma;

			k = (k + 1) % x2;
			t--;
			if (t%x2 == 0 && t != t1 && w < (y1 - 1)) w++;
			//t--;
			if (t == 0) break;
		}

		return Matrix(y1, x2, newArray);
	}
}

Matrix Matrix::operator *(const double& rhs) {
		int y1 = this->liczba_macierzy;
		int x1 = this->liczba_elementow;

		double** newArray = new double*[y1];
		for (int i = 0; i < y1; i++)
			newArray[i] = new double[x1];

		for (int y = 0; y < y1; y++) {
			for (int x = 0; x < x1; x++) {
				newArray[y][x] = this->arrays[y][x] * rhs;
			}
		}

		return Matrix(y1, x1, newArray);
}

Matrix::Matrix(unsigned int liczba_macierzy, unsigned int liczba_elementow, double** arrays) {
	this->liczba_elementow = liczba_elementow;
	this->liczba_macierzy = liczba_macierzy;
	this->arrays = arrays;
}

Matrix::Matrix(unsigned int liczba_macierzy, unsigned int liczba_elementow) {
	this->liczba_macierzy = liczba_macierzy;
	this->liczba_elementow = liczba_elementow;
	this->initMatrix();
}

Matrix::Matrix(unsigned int liczba_elementow) {
	this->liczba_elementow = liczba_elementow;
	this->liczba_macierzy = 0;
	this->arrays = nullptr;
}

Matrix::Matrix() {
	this->liczba_elementow = 0;
	this->liczba_macierzy = 0;
	this->arrays = nullptr;
}

void Matrix::initMatrix() {
	this->arrays = new double*[liczba_macierzy];
	for (int i = 0; i < liczba_macierzy; i++) {
		arrays[i] = new double[liczba_elementow];
	}
}

void Matrix::add(unsigned int count, ...) {
	this->liczba_macierzy = count;
	this->initMatrix();
	double* arg = new double[this->liczba_elementow];
	unsigned int j = 0;

	va_list argument;
	__crt_va_start(argument, count);
	for (int i = 0; i < count; i++) {

		arg = __crt_va_arg(argument, double*);
		for (int i = 0; i < this->liczba_elementow; i++) {
			this->arrays[j][i] = arg[i];
		}
		j++;

	}
	__crt_va_end(argument);
}

void Matrix::add(const double *macierz, ...) {
	unsigned int j = 0;
	va_list arguments;
	for (__crt_va_start(arguments, macierz); macierz != NULL && j < this->liczba_macierzy; macierz = __crt_va_arg(arguments, const double*)) {

		for (int i = 0; i < this->liczba_elementow; i++) {
			this->arrays[j][i] = macierz[i];
		}
		j++;

	} __crt_va_end(arguments);
}

void Matrix::add(string macierz) {
	//"[1,2,1.4][3,9,3][41,2,55]"

	//Podzia³ na pojedyñcze macierze (fragmenty)
	unsigned int len = macierz.length();
	string fragment = "";
	vector <string> fragmenty;

	for (unsigned int i = 0; i < len; i++) {
		if (macierz[i] == '[') {
			fragment = "";
			do {
				i++;
				if (macierz[i] == ']') break;

				fragment += macierz[i];
			} while (i < len - 1);
			fragmenty.push_back(fragment);
		}
	}

	//cout << "Fragmenty: "; for (int i = 0; i < fragmenty.size(); i++) cout << fragmenty[i] << " ; "; cout << endl;

	this->liczba_macierzy = fragmenty.size();

	//Podzia³ fragmentów na pojedyñcze elementy i sprawdzenie iloœci elementów
	vector <double> wartosci;
	short elementy = 0;

	for (int i = 0; i < this->liczba_macierzy; i++) {
		len = fragmenty[i].size(); //d³ugoœæ stringa
		
		for (int j = 0; j < len; j++) {
			fragment = "";
			while (fragmenty[i][j] != ',' && j < len) {
				fragment += fragmenty[i][j];
				j++;
			}
			wartosci.push_back(stod(fragment, string::size_type()));
		}

		if (i == 0) elementy = wartosci.size();
	}

	//cout << "Wartosci: "; for (int i = 0; i < wartosci.size(); i++) cout << wartosci[i] << " ; "; cout << endl;

	this->liczba_elementow = elementy;

	//Inicjowanie nowej macierzy typu Matrix
	this->initMatrix();

	//Przenoszenie elementów z vector do macierzy typu Matrix
	int i = 0;
	for (int y = 0; y < this->liczba_macierzy; y++) {
		for (int x = 0; x < this->liczba_elementow; x++) {
			this->arrays[y][x] = wartosci[i];
			i++;
		}
	}

	//czyszczenie pamiêci
	fragmenty.clear();
	wartosci.clear();
}

Matrix Matrix::print() {
	for (int i = 0; i < this->liczba_macierzy; i++) {
		cout << "[";
		for (int j = 0; j < this->liczba_elementow; j++) {
			cout << " " << arrays[i][j];
		}
		cout << " ]" << endl;
	}
	return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
}

unsigned int power(short inp) {
	unsigned int output = 1;
	for (short i = 0; i < inp; i++)
		output *= 10;
	return output;
}

Matrix Matrix::print(short roundness) {
	int pomocnicza = 0;
	roundness = (roundness < 5 ? roundness : 4);
	for (int i = 0; i < this->liczba_macierzy; i++) {
		cout << "[";
		for (int j = 0; j < this->liczba_elementow; j++) {
			if (roundness != 0) {
				pomocnicza = (float)arrays[i][j] * power(roundness);
				cout << " " << ((float)pomocnicza / power(roundness));
			}
			else
				cout << " " << round(arrays[i][j]);
		}
		cout << " ]" << endl;
	}
	return Matrix(this->liczba_macierzy, this->liczba_elementow, this->arrays);
}



class NeuralNetwork{
public:
	friend class NeuNet;
	NeuralNetwork(unsigned int,unsigned int, int seed = 1);
	NeuralNetwork(unsigned int neuron_inputs = 0);
	double drand(double, double);
	void train(Matrix, Matrix, unsigned int);
	Matrix think(Matrix);
	void print_synaptic_weights();
	void print_names();
	void add_names(string);
	void print_classified();

	//eksperymentalne
	Matrix dot(Matrix,Matrix);
protected:
	unsigned int neuron_count;
	unsigned int neuron_inputs;
	Matrix* synaptic_weights;
	Matrix wynik;
	vector <string> nazwy;
};

Matrix NeuralNetwork::dot(Matrix lhs, Matrix rhs) {
	if (!(lhs.liczba_elementow == rhs.liczba_macierzy)) {
		cout << "Nie mozna pomnozyc tych macierzy" << endl;
		return Matrix(lhs.liczba_macierzy, lhs.liczba_elementow, lhs.arrays);
	}
	else {
		int y1 = lhs.liczba_macierzy;
		int x1 = lhs.liczba_elementow;
		int y2 = rhs.liczba_macierzy;
		int x2 = rhs.liczba_elementow;

		double** newArray = new double*[y1];
		for (int i = 0; i < y1; i++) {
			newArray[i] = new double[x2];
		}

		double suma = 0;
		int w, k, t;
		w = k = 0;
		t = rhs.liczba_elementow * lhs.liczba_macierzy;
		int t1 = t;

		while (t > 0) {
			suma = 0;
			for (int i = 0; i < x1; i++) {
				suma += lhs.arrays[w][i] * rhs.arrays[i][k];
			}

			newArray[w][k] = suma;

			k = (k + 1) % x2;
			t--;
			if (t%x2 == 0 && t != t1 && w < (y1 - 1)) w++;
			//t--;
			if (t == 0) break;
		}

		//czyszczenie
		for (int i = 0; i < y2; i++) delete[] rhs.arrays[i];
		for (int i = 0; i < y1; i++) delete[] lhs.arrays[i];
		delete[] lhs.arrays;
		delete[] rhs.arrays;


		return Matrix(y1, x2, newArray);
	}
}

void NeuralNetwork::print_classified() {
	unsigned int klasy = this->neuron_count;

	for (int i = 0; i < klasy; i++) {
		if (this->wynik.getArray()[0][i] >= 0.5) cout << "[ " << this->nazwy[i] << " ] ";
	}
	cout << endl;
}

void NeuralNetwork::print_names() {
	cout << "Stored names: ";
	for (int i = 0; i < this->nazwy.size(); i++) {
		cout << "[ " << this->nazwy[i] << " ] ";
	}
	cout << endl;
}

void NeuralNetwork::add_names(string input) {
	//Podzia³ na pojedyñcze macierze (fragmenty)
	unsigned int len = input.length();
	string fragment = "";

	for (unsigned int i = 0; i < len; i++) {
		if (input[i] == '[') {
			fragment = "";
			do {
				i++;
				if (input[i] == ']') break;

				fragment += input[i];
			} while (i < len - 1);
			this->nazwy.push_back(fragment);
		}
	}
}

void NeuralNetwork::train(Matrix training_inputs, Matrix training_outputs, unsigned int iterations) {
	Matrix output(3);
	Matrix error(3);
	Matrix adjustment(3);

	unsigned int modulo = 5 * (iterations / 100);

	cout << " [ ";

	for (unsigned int i = 0; i < iterations; i++) {

		if (i%modulo == 0) cout << (i * 100) / iterations << "%  ";

		output = this->think(training_inputs);

		error = training_outputs - output;

		//adjustment = training_inputs.T() * (output.sigmoid_derivative() *= error);
		adjustment = dot(training_inputs.T(), (output.sigmoid_derivative() *= error));  //zjada mniej ramu


		synaptic_weights += adjustment;
	}

	cout << " 100% ] " << endl;
}

Matrix NeuralNetwork::think(Matrix inputs) {
	this->wynik = (inputs * synaptic_weights).sigmoid();
	return this->wynik;
}

NeuralNetwork::NeuralNetwork(unsigned int neuron_inputs, unsigned int neuron_count , int seed) {
	srand(seed);
	this->neuron_count = neuron_count;
	this->neuron_inputs = neuron_inputs;

	this->synaptic_weights = new Matrix(neuron_inputs);

	string weights = "";
	string liczba = "";
	for (int j = 0; j < neuron_count; j++) {
		weights += '[';
		for (int i = 0; i < neuron_inputs; i++) {
			liczba = to_string(this->drand(0, 2) - 1);
			for (int z = 0; z < liczba.size(); z++)
				weights += liczba[z];
			weights += (i == neuron_count ? ' ' : ',');
		}
		weights += ']';
	}

	//cout << "Wagi: " << weights << endl;

	this->synaptic_weights->add(weights);

	synaptic_weights = synaptic_weights->t();
}

NeuralNetwork::NeuralNetwork(unsigned int neuron_inputs) {
	unsigned int neuron_count = 1;
	int seed = 1;
	srand(seed);
	this->neuron_count = neuron_count;
	this->neuron_inputs = neuron_inputs;
	this->synaptic_weights = new Matrix(neuron_inputs);

	double* arr = new double[neuron_inputs];
	for (int i = 0; i < neuron_inputs; i++) {
		arr[i] = this->drand(0, 2) - 1;
	}
	this->synaptic_weights->add(neuron_count, arr);

	delete[] arr;

	synaptic_weights = synaptic_weights->t();
}

double NeuralNetwork::drand(double min, double max) {
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

void NeuralNetwork::print_synaptic_weights() {
	synaptic_weights->print();
}

class NeuNet {
public:
	//NeuNet(string, string, string, unsigned int, bool);
	NeuNet(bool sigm = 1);
	void input(string);
	void output(string);
	void labels(string);
	void iterations(unsigned int);
	bool Setup();
	bool Train();
	void Think(string);
protected:
	NeuralNetwork neural_net;
	Matrix training_inputs;
	Matrix training_outputs;
	bool sigm;
	bool setup;
	unsigned int iteration;
	string names;
	string ID;
};

void NeuNet::Think(string data) {
	cout << ID << "Considering new situation: " << data << endl;
	Matrix nowa;
	nowa.add(data);
	if (this->sigm) nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	neural_net.print_classified();
	cout << endl;
}

bool NeuNet::Train() {
	cout << ID << "Starting Training..." << endl;
	if (setup) {
		std::clock_t start;
		double durationTh;
		cout << ID << "Iterations: " << this->iteration << endl;
		start = std::clock();

		neural_net.train(training_inputs, training_outputs, this->iteration);

		durationTh = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		cout << ID << "Succeeded in time: [" << durationTh << "] s" << endl << endl;

		cout << ID << "New synaptic weights after training: " << endl;
		neural_net.print_synaptic_weights();
		cout << endl << endl;
	}
	else {
		cout << ID << "Training failed. (Try Setup first)" << endl;
		return 0;
	}
}

bool NeuNet::Setup() {
	cout << ID << "Starting Setup..." << endl;
	if (!(this->training_inputs.arrays == nullptr || this->training_outputs.arrays == nullptr || this->iteration == 0 || this->neural_net.nazwy.size() == 0)) {
		cout << ID << "Setting up NeuralNetwork" << endl;
		NeuralNetwork neur(this->training_inputs.liczba_elementow, this->training_outputs.liczba_elementow);
		neur.add_names(this->names);
		this->neural_net = neur;

		cout << ID << "Random starting synaptic weights:" << endl;
		this->neural_net.print_synaptic_weights();
		cout << endl;
		setup = 1;
		return 1;
	}
	else {
		cout << ID << "Setup failed" << endl;
		return 0;
	}
}

NeuNet::NeuNet(bool sigm) {
	this->ID = "[NeuNet] ";
	cout << ID << "Created instance" << endl;
	this->sigm = sigm;
	iteration = 0;
	setup = 0;
}

void NeuNet::input(string training_data) {
	cout << ID << "Adding training input" << endl;
	this->training_inputs.add(training_data);

	if (this->sigm) {
		cout << ID << "Using sigmoid function on training input" << endl;
		this->training_inputs = this->training_inputs.sigmoid();
	}
}

void NeuNet::output(string training_output) {
	cout << ID << "Adding training output" << endl;
	this->training_outputs.add(training_output);
	this->training_outputs = this->training_outputs.T();
}

void NeuNet::labels(string names) {
	cout << ID << "Adding labels" << endl;
	this->neural_net.add_names(names);
	this->names = names;
}

void NeuNet::iterations(unsigned int iter) {
	this->iteration = iter;
}




void Test1() {
	NeuNet n;

	n.input("[1,5][1,6][1,4.5][1,5.25][1,4] [0,12][0,0.14][0,13][0,13.5][0,17]");
	n.output("[1,1,1,1,1,0,0,0,0,0] [0,0,0,0,0,1,1,1,1,1]");
	n.labels("[japko][pomaranicz]");
	n.iterations(1000000);

	n.Setup();
	n.Train();

	n.Think("[1,3]");
	n.Think("[0,5]");
	n.Think("[0,12]");
	n.Think("[1,32]");
	n.Think("[0,1]");
}

void Test2() {
	NeuNet n;

	n.input("[0,0,1,0][0,0,0,1] [1,1,1,0][0,1,1,1] [1,1,1,1] [0,0,1,1][0,1,1,0][0,1,0,0][1,0,1,1][1,1,0,1]");
	n.output("[1,1,1,1,0,0,0,1,1,1] [0,0,1,1,1,1,1,0,1,1] [0,0,0,0,1,0,0,0,0,0]");
	n.labels("[jeden][dwa][cztery]");
	n.iterations(100000);

	n.Setup();
	n.Train();

	n.Think("[1,0,1,1]");
	n.Think("[1,1,0,0]");
	n.Think("[0,1,0,0]");
	n.Think("[1,0,0,0]");
	n.Think("[1,1,1,1]");
	n.Think("[0,0,1,0]");
}

void Test3() {
	NeuNet n;
	n.input("[3,10][2.4,12][3.1,11] [9,0.1][9.5,3][9.2,2]");
	n.output("[1,1,1,0,0,0] [0,0,0,1,1,1]");
	n.labels("[japko][pomaranicz]");
	n.iterations(1000000);

	n.Setup();
	n.Train();

	n.Think("[1,13]");
	n.Think("[5,15]");
	n.Think("[10,1]");
}

int main() {
	cerr.sync_with_stdio(false);

	//Test1();  //zu¿ywa najwiêcej ramu
	//Test2();

	Test3();

	cout << endl;
	_getch();
}