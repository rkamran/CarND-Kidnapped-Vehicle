/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 20;

	for(int i=0; i<num_particles; i++){
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
		weights.push_back(p.weight);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);



	for(int i=0; i<num_particles; i++){
		Particle p = particles[i];

		double sin_theta = sin(p.theta);
		double cos_theta = cos(p.theta);
		double sin_theta_plus_yaw = sin(p.theta + yaw_rate*delta_t);
		double cos_theta_plus_yaw = cos(p.theta + yaw_rate*delta_t);

		if (fabs(yaw_rate) < 0.0001){
		    p.x = dist_x(gen) + velocity*delta_t*cos_theta;
		    p.y = dist_y(gen) + velocity*delta_t*sin_theta;
		}
		else{
			p.x = dist_x(gen) + p.x + ((velocity/yaw_rate)*(sin_theta_plus_yaw - sin_theta));
			p.y = dist_y(gen) + p.y + ((velocity/yaw_rate)*(cos_theta - cos_theta_plus_yaw));
			p.theta = dist_theta(gen) + p.theta + yaw_rate*delta_t;
		}
		particles[i] = p;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i=0; i<predicted.size(); ++i){

		double distance = 0;
		int nearest_index = -1;

		for (int j=0; j< observations.size(); ++j){
			if(j == 0){
				distance = sqrt(pow((predicted[i].x - observations[j].x), 2) + pow((predicted[i].y - observations[j].y), 2));
			}
			else{
				double new_distance = sqrt(pow((predicted[i].x - observations[j].x), 2) + pow((predicted[i].y - observations[j].y), 2));
				if(new_distance < distance){
					distance = new_distance;
					nearest_index = j;
				}
			}
		}
		if (nearest_index > -1){
			predicted[i].x = observations[nearest_index].x;
			predicted[i].y = observations[nearest_index].y;
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//-> Calcualte gaussian
	auto gaussian =  [std_landmark] (double mu_x, double mu_y, double x_obs, double y_obs){
		double gauss_norm= 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);
		double exponent= ((x_obs - mu_x) * (x_obs - mu_x)) / (2 * (std_landmark[0] * std_landmark[0])) +
						 ((y_obs - mu_y) * (y_obs - mu_y)) / (2 * (std_landmark[1] * std_landmark[1]));

		return gauss_norm * exp(-exponent);
	};

	//-> Transform observations for a particle to Map Cooridnate
	auto transform = [observations] (Particle p){
		vector<LandmarkObs> transformed_obs;

		for (LandmarkObs obs : observations){
			LandmarkObs tr_obs;
			tr_obs.x = p.x + (cos(p.theta) * obs.x) + (sin(p.theta) * obs.y);
			tr_obs.y = p.y - (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
			tr_obs.id = p.id;
			transformed_obs.push_back(tr_obs);
		}
		return transformed_obs;
	};

	//-> Associate
	auto associate = [this, map_landmarks] (vector<LandmarkObs> tr_obs, Particle& p){

		std::vector<LandmarkObs> predictions;
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		for(int i=0; i< tr_obs.size(); ++i){
			double min_distance = std::numeric_limits<float>::max();
			int nearest_id = -1;
			double nearest_x = 0.0;
			double nearest_y = 0.0;

			for (int j=0; j< map_landmarks.landmark_list.size(); ++j){
				double this_distance = dist(tr_obs[i].x, tr_obs[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);

				if (this_distance < min_distance){
					min_distance = this_distance;

					nearest_x = map_landmarks.landmark_list[j].x_f;
					nearest_y = map_landmarks.landmark_list[j].y_f;
					nearest_id = map_landmarks.landmark_list[j].id_i;
				}
			}

			sense_x.push_back(nearest_x);
			sense_y.push_back(nearest_y);
			associations.push_back(nearest_id);
		}
		p.sense_x.clear(); p.sense_x = sense_x;
		p.sense_y.clear(); p.sense_y = sense_y;
		p.associations.clear(); p.associations = associations;
	};

	weights.clear();

	//For each particle
	for(int i=0; i<num_particles; i++){
		//Transform observations to map cordinate
		vector<LandmarkObs> tr_obs = transform(particles[i]);

		//Assoctiate
		associate(tr_obs, particles[i]);

		double prob = 1.0;

		for (int ii=0; ii < particles[i].associations.size(); ++ii){
			double mark_x = particles[i].sense_x[ii];
			double mark_y = particles[i].sense_y[ii];
			double weight = gaussian(mark_x, mark_x, tr_obs[ii].x, tr_obs[ii].y);
			if(weight > 0)
				prob *= weight;
		}
		particles[i].weight = prob;
		weights.push_back(prob);
		cout << "Weight 0 -- "<< weights[i] <<endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	auto random_gen = [] (){
		return ((double) rand() / (RAND_MAX)) + 1;
	};

	vector<Particle> resampled;
	int index = rand() % num_particles;
	double beta = 0.0;
	double mw = *max_element(weights.begin(), weights.end());

	for(int i=0; i<num_particles; i++){
		beta+= random_gen() * 2.0 * mw;
		while (beta > weights[index]){
			beta -= weights[index];
			index = (index +1) % num_particles;
		}
		resampled.push_back(particles[index]);
	}
	particles = resampled;
	weights.clear();
	for(int i=0; i<num_particles; i++)
		weights.push_back(particles[i].weight);
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
