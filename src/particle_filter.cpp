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

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 100;

	for(int i=0; i<num_particles; i++){
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(1.0);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);


	for(int i=0; i<num_particles; i++){
		if (fabs(yaw_rate) < 0.001){
			particles[i].x += (velocity * delta_t * cos(particles[i].theta)) + dist_x(gen);
			particles[i].y += (velocity * delta_t * sin(particles[i].theta)) + dist_y(gen);
			particles[i].theta += (yaw_rate * delta_t) + dist_theta(gen);
		}
		else{
			particles[i].x +=  (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y +=  (velocity/yaw_rate) * (cos(particles[i].theta ) - cos(particles[i].theta + yaw_rate)) + dist_y(gen);
			particles[i].theta += (yaw_rate * delta_t) + dist_theta(gen);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i=0; i< observations.size(); i++){

		double min_distance = numeric_limits<double>::max();
		int pred_id = -1;

		LandmarkObs obs = observations[i];

		for (int j=0; j<predicted.size(); j++){

			LandmarkObs pred = predicted[j];

			double this_distance = dist(obs.x, obs.y, pred.x, pred.y);

			if(this_distance < min_distance){
				min_distance = this_distance;
				pred_id = j;
			}
		}
		observations[i].id = pred_id;
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


	for(int i=0; i<num_particles; i++){

		//-->1. All map landmark within in range
		vector<LandmarkObs> valid_landmarks;

		//-->2. This particle
		Particle p = particles[i];

		//-->3. Filter out far away map landmarks
		cout << "Before " << valid_landmarks.size();
		for(Map::single_landmark_s lm : map_landmarks.landmark_list){
			if(fabs(lm.x_f - p.x) > sensor_range || fabs(lm.y_f - p.y) > sensor_range)
				continue;
			valid_landmarks.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
		}
		cout << "After " << valid_landmarks.size();

		//-->4. Transform observations to map co-ordinate for this particle
		vector<LandmarkObs> transformed_obs;

		for (LandmarkObs obs : observations){
			LandmarkObs tr_obs;
			tr_obs.x = p.x + obs.x * cos(p.theta) - obs.y * sin(p.theta);
			tr_obs.y = p.y + obs.x * sin(p.theta) + obs.y * cos(p.theta);
			tr_obs.id = obs.id;
			transformed_obs.push_back(tr_obs);
		}

		//-->5. Associate valid map landmarks to observations
		dataAssociation(valid_landmarks, transformed_obs);

		//-->6. Calcualte particle weight
		double prob = 1.0;

		for(int j=0; j<transformed_obs.size(); j++){
			double mu_x, mu_y;
			double x_obs, y_obs;

			x_obs = transformed_obs[j].x;
			y_obs = transformed_obs[j].y;

			if(transformed_obs[j].id != -1){
				LandmarkObs lm = valid_landmarks[transformed_obs[j].id];
				mu_x = lm.x;
				mu_y = lm.y;
			}

			double gauss_norm= 2.0 * M_PI * std_landmark[0] * std_landmark[1];
			double exponent= (pow((x_obs - mu_x),2)) / (2 * pow(std_landmark[0],2)) +
							(pow((y_obs - mu_y),2)) / (2 * pow(std_landmark[1],2));
			prob *=  exp(-exponent) / gauss_norm;
		}
		particles[i].weight = prob;
		weights[i] = prob;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::discrete_distribution<int> d(weights.begin(), weights.end());
	std::vector<Particle> new_particles;

	for(unsigned i = 0; i < num_particles; i++)
	{
		auto ind = d(gen);
		new_particles.push_back(std::move(particles[ind]));
	}
	particles = std::move(new_particles);
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
