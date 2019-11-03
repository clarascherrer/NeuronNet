#include "network.h"
#include "random.h"

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}


std::pair<size_t, double> Network::degree(const size_t& n) const
{
	std::pair<size_t, double> deg;
	size_t number_connections(neighbors(n).size());
	double tot_intensity;
	
	for(auto connexion : neighbors(n))
	{
		tot_intensity += connexion.second;
	}
	
	deg = std::make_pair(number_connections, tot_intensity);
	
	return deg;
}

//std::set is a container that stores unique elements following a specific order : 
std::set<size_t> Network::step(const std::vector<double>& thalamic_inputs)
{
	std::set<size_t> firing_neurons;
	
	//Fills the set with the indices of firing neurons
	for(size_t n(0); n < neurons.size(); ++n) {
		if(neurons[n].firing()) {
			firing_neurons.insert(n);
			neurons[n].reset(); //after a neuron fired it is reset
		
		} else { //only non firing neurons are reset
		
		//Evolution of the neuron's intensity :
		
		std::vector<std::pair<size_t, double>> connected_neurons(neighbors(n));
		double noise(thalamic_inputs[n]);
		double new_intensity;
		
		//each neuron receives an input from firing connected neurons :
		
		for(auto neighbor : connected_neurons) {
			if(!neurons[neighbor.first].is_inhibitory()) {
				new_intensity += 0.5*neighbor.second;
			}
			
			else if(neurons[neighbor.first].is_inhibitory()) {
				new_intensity -= neighbor.second;
			}
		}
		
		if(neurons[n].is_inhibitory()) {
			noise *= 0.4;
		}
		
	//Intakes the new calculated input in neuron n
	neurons[n].input(noise + new_intensity);
	
	//Update of the neuron's potential and recovery according to Izhikevich equations :
	neurons[n].step();
		}
	}
	
	return firing_neurons;
}


//Finds the list of neurons connected to the neuron studied of index n
std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& n) const
{
	std::vector<std::pair<size_t, double> > connected_neurons;
	
	//For 2 neurones to be connected, the index of the first neuron (n) must be equal to
	// the first value in the pair of the linkmap corresponding to the same index.
	//Using the map iterators optimises the program, making it much faster than iterating on the linkmap
	
	std::map<std::pair<size_t, size_t>, double>::const_iterator low_iterator;
	std::map<std::pair<size_t, size_t>, double>::const_iterator up_iterator;
	
	low_iterator = links.lower_bound(std::make_pair(n,0)); //pair between the neuron n and any other neuron : (n,0), the second index doesn't matter
	up_iterator = links.upper_bound(std::make_pair(n,0));
	
	for(low_iterator; low_iterator != up_iterator; ++low_iterator) {
		connected_neurons.push_back(std::make_pair(low_iterator->first.second, low_iterator->second)); //first.second = index of neuron connected to neuron n and low_iterator->second = intensity of connection 
	}
	
	return connected_neurons;
}
