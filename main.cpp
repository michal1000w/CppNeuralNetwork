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
	friend Matrix* operator += (const Matrix* lhs, const Matrix& rhs);
	Matrix(unsigned int, unsigned int, double**);
	Matrix(unsigned int, unsigned int);
	Matrix(unsigned int);
	Matrix();
	void add(const double*, ...);
	void add(unsigned int, ...);
	void add(string);
	void print();
	void print(short);
	Matrix T();
	Matrix* t();

	//funkcjonalno�� pod AI
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
private:
	void initMatrix();
	unsigned int liczba_macierzy;
	unsigned int liczba_elementow;
	double** arrays;
};

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

	//Podzia� na pojedy�cze macierze (fragmenty)
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

	//Podzia� fragment�w na pojedy�cze elementy i sprawdzenie ilo�ci element�w
	vector <double> wartosci;
	short elementy = 0;

	for (int i = 0; i < this->liczba_macierzy; i++) {
		len = fragmenty[i].size(); //d�ugo�� stringa
		
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

	//Przenoszenie element�w z vector do macierzy typu Matrix
	int i = 0;
	for (int y = 0; y < this->liczba_macierzy; y++) {
		for (int x = 0; x < this->liczba_elementow; x++) {
			this->arrays[y][x] = wartosci[i];
			i++;
		}
	}

	//czyszczenie pami�ci
	fragmenty.clear();
	wartosci.clear();
}

void Matrix::print() {
	for (int i = 0; i < this->liczba_macierzy; i++) {
		cout << "[";
		for (int j = 0; j < this->liczba_elementow; j++) {
			cout << " " << arrays[i][j];
		}
		cout << " ]" << endl;
	}
}

unsigned int power(short inp) {
	unsigned int output = 1;
	for (short i = 0; i < inp; i++)
		output *= 10;
	return output;
}

void Matrix::print(short roundness) {
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
}



class NeuralNetwork{
public:
	NeuralNetwork(unsigned int,unsigned int, int seed = 1);
	NeuralNetwork(unsigned int);
	double drand(double, double);
	void train(Matrix, Matrix, unsigned int);
	Matrix think(Matrix);
	void print_synaptic_weights();
protected:
	unsigned int neuron_count;
	unsigned int neuron_inputs;
	Matrix* synaptic_weights;
};

void NeuralNetwork::train(Matrix training_inputs, Matrix training_outputs, unsigned int iterations) {
	Matrix output(3);
	Matrix error(3);
	Matrix adjustment(3);
	for (unsigned int i = 0; i < iterations; i++) {
		output = this->think(training_inputs);

		//cout << "Output: " << endl;
		//output.print();

		error = training_outputs - output;

		//cout << "Error: " << endl;
		//error.print();


		adjustment = training_inputs.T() * (output.sigmoid_derivative() *= error);

		//cout << "Adjustment" << endl;
		//adjustment.print();

		synaptic_weights += adjustment;
	}
}

Matrix NeuralNetwork::think(Matrix inputs) {
	return (inputs * synaptic_weights).sigmoid();
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

	synaptic_weights = synaptic_weights->t();
}

double NeuralNetwork::drand(double min, double max) {
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

void NeuralNetwork::print_synaptic_weights() {
	synaptic_weights->print();
}


int main() {
#define a new double
	NeuralNetwork neural_net(2,2);

	cout << "Random starting synaptic weights: " << endl;
	neural_net.print_synaptic_weights();
	cout << endl;


	Matrix training_inputs;
	training_inputs.add("[1,5][1,6][1,4.5][1,5.25][1,4] [0,12][0,0.14][0,13][0,13.5][0,17]");
	training_inputs = training_inputs.sigmoid(); ///test
	Matrix training_outputs;
	training_outputs.add("[1,1,1,1,1,0,0,0,0,0] [0,0,0,0,0,1,1,1,1,1]");
	training_outputs = training_outputs.T();

	//////////Trening
	std::clock_t start;
	double durationTh;
	start = std::clock();
	cout << "Rozpoczynam trening..." << endl;

		neural_net.train(training_inputs, training_outputs, 100000);

	durationTh = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "Zakonczono pomyslnie w czasie: [" << durationTh << "] s" << endl << endl;
	//////////Koniec

	cout << "New synaptic weights after training: " << endl;
	neural_net.print_synaptic_weights();
	cout << endl;

	cout << "Considering new situation [1,3]" << endl;
	Matrix nowa;
	nowa.add("[1,3]");
	nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	cout << endl;

	cout << "Considering [0,5]" << endl;
	nowa.add("[0,5]");
	nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	cout << endl;

	cout << "Considering [0,12]" << endl;
	nowa.add("[0,12]");
	nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	cout << endl;

	cout << "Considering [1,32]" << endl;
	nowa.add("[1,32]");
	nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	cout << endl;

	cout << "Considering [0,1]" << endl;
	nowa.add("[0,1]");
	nowa = nowa.sigmoid();
	neural_net.think(nowa).print(0);
	cout << endl;


	_getch();
}