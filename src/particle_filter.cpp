/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles  
  std::default_random_engine gen;

  // Create normal distributions for y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta;   
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);   
    
	Particle p_temp;
    p_temp.id = i;
    p_temp.x = sample_x;
    p_temp.y = sample_y;
    p_temp.theta = sample_theta;
    p_temp.weight = 1.0;
    particles.push_back(p_temp);
    
    weights.push_back(1.0); 
    // Print your samples to the terminal.
    std::cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " 
              << sample_theta << std::endl;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;

  // Create normal distributions for y and theta
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta, x, y, theta;   
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);   
      
    x = particles[i].x;
    y = particles[i].y;
    theta = particles[i].theta;
    
    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x = x + velocity*delta_t*cos(theta) + sample_x;
      particles[i].y = y + velocity*delta_t*sin(theta) + sample_y;
    } else {
      particles[i].x = x + (velocity/yaw_rate)*(sin(theta + yaw_rate*delta_t) - sin(theta)) + sample_x;
      particles[i].y = y + (velocity/yaw_rate)*(cos(theta) - cos(theta + yaw_rate*delta_t)) + sample_y;
      particles[i].theta = theta + yaw_rate*delta_t + sample_theta;
    }        
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for (unsigned int i=0; i<observations.size(); i++) {
  	double min = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
  	int associate_id = predicted[0].id;
    for (unsigned int j=1; j<predicted.size(); j++) {
      double current = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (current < min) {
       	min = current;
        associate_id = predicted[j].id;
      }
    }
    observations[i].id = associate_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double weights_sum = 0.0;
  std::vector<LandmarkObs> predictions;
  std::vector<LandmarkObs> observ_map;
  
  for (int i = 0; i < num_particles; ++i) {
  	particles[i].associations.clear();
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
  }
  
  for (int i=0; i<num_particles; i++) {
  	
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++) {
      LandmarkObs temp;
      temp.x = map_landmarks.landmark_list[j].x_f;
      temp.y = map_landmarks.landmark_list[j].y_f;
      temp.id = map_landmarks.landmark_list[j].id_i; 
      if ((fabs(temp.x-particles[i].x)<=sensor_range) && (fabs(temp.y-particles[i].y)<=sensor_range)) {
        predictions.push_back(temp);
      }
    }
  
    for (unsigned int j = 0; j < observations.size(); j++) {
      LandmarkObs temp;
      temp.x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
      temp.y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
      //temp.id = observations[j].id; 
      observ_map.push_back(temp);
    }
  
  	dataAssociation(predictions, observ_map);
    
    for (unsigned int j = 0; j < observ_map.size(); j++) {
      particles[i].associations.push_back(observ_map[j].id);
      particles[i].sense_x.push_back(observ_map[j].x);
      particles[i].sense_y.push_back(observ_map[j].y);
      std::cout << "ID " << j + 1 << ": " << observ_map[j].id << "\t";
    }
    
    
    particles[i].weight = 1.0;
    for (unsigned int j = 0; j < observ_map.size(); j++) {
      double observ_x, observ_y, predict_x, predict_y;
      observ_x = observ_map[j].x;
      observ_y = observ_map[j].y;
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (observ_map[j].id == predictions[k].id) {
          predict_x = predictions[k].x;
          predict_y = predictions[k].y;
         // break;
        }
      }
      double observ_w = multiv_prob(std_landmark[0], std_landmark[1], observ_x, observ_y, predict_x, predict_y);
      particles[i].weight *= observ_w;
    }
    
    weights[i] = particles[i].weight;
    weights_sum += weights[i];
    predictions.clear();
    observ_map.clear();
    std::cout << std::endl << predictions.size() << std::endl;
  }
  
  for (int i=0; i<num_particles; i++) {
    particles[i].weight /= weights_sum;
    weights[i] /= weights_sum;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<int> distribution (weights.begin(), weights.end());
  std::vector<Particle> sample_particles;
  for (int i = 0; i < num_particles; ++i) {
  	sample_particles.push_back(particles[distribution(gen)]);
  }
  particles = sample_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}