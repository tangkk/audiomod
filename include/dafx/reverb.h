/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

# pragma once

#include "modbase.h"

#include <mutex>

class revmodel;

class reverb : public modbase
{   
public:
    reverb(int sampleRate, int numChannels, float roomsize = 0.8f, float damp = 0.9f, float width = 0.5f, float dry = 1.f, float wet = 0.2f);
    ~reverb();

    /**
     * set mod params
     * @param params the key-val params to be set
     */ 
    void setParams(std::map<std::string, float> params) {

    }

    /**
     * get mod params
     * @param params the key-val params to be returned
     */ 
    void getParams(std::map<std::string, float> &params) {

    }

    void	setroomsize(float value);
    float	getroomsize() const;
    void	setdamp(float value);
    float	getdamp() const;
    void	setwet(float value);
    float	getwet() const;
    void	setdry(float value);
    float	getdry() const;
    void	setwidth(float value);
    float	getwidth() const;
    void	setmode(float value);
    float	getmode() const;


    void processBlock (float *const * bufferData, int num_samples);
private:
    revmodel *rev;
    std::mutex coffupdate_mutex;

};
