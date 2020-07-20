
Code and data for

#### Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, 8:e46331. ####

Behavioral data and model fits are available at https://doi.org/10.6084/m9.figshare.7268558 under a CC-BY 4.0 license.

To fit the models, install the HDDM package (http://ski.clps.brown.edu/hddm_docs/index.html). Then run the HDDM models using b1_HDDM_run.py (the models are specified in hddm_models.py). Easiest is to use a batch job submission system, and do e.g.
<code>
python b1_HDDM_run.py -r 1 -d $d -v $v -i $i -s $s
</code>
where -d = 0-6 (datasets), -v = 0-11 (versions of the model), -i = 0-30 (traces, can be changed to whatever the number of cores on a node) and -i = 5.000, the number of samples per trace.

To reproduce all main figures, see <code>plot_all.m</code> which has the overview of the scripts that are called to generate each figure. 

The <code>extended_models</code> folder contains Matlab code to fit the models in Figure 6.

The <code>simulations</code> folder contains Python code to generate Supplementary Figure 8.

For questions, @AnneEUrai / anne.urai@gmail.com.


#### LICENSE
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. _If you use the Software for your own research, please cite the paper._

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
