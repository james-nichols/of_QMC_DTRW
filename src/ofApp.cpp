#include "ofApp.h"
#include <math.h>
#include <time.h>

//--------------------------------------------------------------
void ofApp::setup() {
     
    bRun = false;

    alpha = 0.8;
    r = 0.5;
    D_alpha = 0.1; 

    T = _T;
    L = _L;
    dX = _DX;
    dT = pow(dX * dX / (2.0 * D_alpha), 1.0 / alpha);
    
    nx = L / dX;
    nt = T / dT;

    nMC = 500000;
    
    ofLog() << "Solving on " << nx << " spatial by " << nt << " time points, dt = " << dT << ", " << nMC << " MC points";
    
    dtrw_solver = new dtrwThread(nx, nt, dX, dT, alpha, r, false);
    mc_solver = new mcThread(nx, nt, dX, dT, nMC, alpha, r, false);
    //qmc_solver = new qmcThread(nx, nt, dX, dT, nMC, alpha, r, false);

    dtrw_solver->start();
    mc_solver->start();
    //qmc_solver->start();

    restart.addListener(this, &ofApp::restartSims);

	gui.setup("MC parameters"); // most of the time you don't need a name
	gui.add(T_slider.setup("End Time", T, 0.01, 5.0 ));
	gui.add(alpha_slider.setup("Alpha", alpha, 0.5, 1.0 ));
	gui.add(D_alpha_slider.setup("D_alpha", D_alpha, 0.01, 1.0 ));
	gui.add(num_mc.setup("Number MC runs", nMC, 1000, 10000000));
	gui.add(restart.setup("Run sim!"));
    gui.setPosition(10, ofGetHeight() - gui.getHeight() - 10);    

    ofEnableSmoothing();
}

//--------------------------------------------------------------
void ofApp::restartSims() {

    T = T_slider;
    alpha = alpha_slider;
    D_alpha = D_alpha_slider;
    dT = pow(dX * dX / (2.0 * D_alpha), 1.0 / alpha);
    
    nt = T / dT;

    nMC = num_mc;

    delete dtrw_solver;
    delete mc_solver;
    //delete qmc_solver;

    ofLog() << "Solving on " << nx << " spatial by " << nt << " time points, dt = " << dT << ", " << nMC << " MC points";
    
    dtrw_solver = new dtrwThread(nx, nt, dX, dT, alpha, r, false);
    mc_solver = new mcThread(nx, nt, dX, dT, nMC, alpha, r, false);
    //qmc_solver = new qmcThread(nx, nt, dX, dT, nMC, alpha, r, false);
    
    dtrw_solver->start();
    mc_solver->start();
    //qmc_solver->start();
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackgroundGradient(ofColor(50), ofColor(0));
   
    dtrw_solver->draw();
    mc_solver->draw();
    //qmc_solver->draw();

    gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
    /*if (key=='r') { 
        bRun = !bRun;
        
        if (bRun) {
            dtrw_solver->start();
            mc_solver->start();
        }
        else {
            dtrw_solver->stop();
            mc_solver->stop();
        }
    }*/
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}
//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}

void dtrwThread::start() 
{ 
    // Mutex blocking is set to true by default
    // It is rare that one would want to use startThread(false).
    startThread();
}

void dtrwThread::stop() 
{
    stopThread();    
}

void dtrwThread::threadedFunction()
{
    while (isThreadRunning()) {
        if (lock()) {
            if (n < nt-1) {
                double out_flux[nx];
            
                n++;
                for (int i=0; i<nx; i++) {
                
                    out_flux[i] = 0.0;
                    
                    for (int m=1; m<=n; m++) {
                        out_flux[i] += U[(n-m)*nx + i] * K[m];
                    }
                    
                    U[n*nx + i] += U[(n-1)*nx + i] - r * out_flux[i];
                    if (i < nx-1)
                        U[n*nx + (i+1)] += r * 0.5 * out_flux[i];
                    else
                        U[n*nx + i] += 0.5 * r * out_flux[i];
                    if (i > 0)
                        U[n*nx + (i-1)] += r * 0.5 * out_flux[i];
                    else
                        U[n*nx + i] += 0.5 * r * out_flux[i];
                }
            }
            else
                stop();

            unlock();
        }
        else {
            // If we reach this else statement, it means that we could not
            // lock our mutex, and so we do not need to call unlock().
            // Calling unlock without locking will lead to problems.
            ofLogWarning("threadedFunction()") << "Unable to lock DTRW mutex.";
        }
    }
}

void dtrwThread::draw() 
{
    if (lock()) {
        ofMesh graph;
        ofMesh graph_double;
        graph.setMode(OF_PRIMITIVE_LINE_STRIP);
        graph_double.setMode(OF_PRIMITIVE_LINE_STRIP);

        double sum = 0.0;
        for (int i=0; i<nx; i++) {
            float x;
            if (centred)
                x = float((ofGetWidth()-2*W_MARGIN) * i) / float(nx) + 0.5/float(nx) + float(W_MARGIN);
            else
                x = float((ofGetWidth()-2*W_MARGIN) * i) / float(nx) + float(W_MARGIN);
            float y = ofGetHeight() - ((ofGetHeight()-2*H_MARGIN) * SCALE * U[n*nx + i] + float(H_MARGIN));
            
            graph.addColor(ofColor::red);
            graph_double.addColor(ofColor::red);
    
            graph.addVertex(ofVec2f(x,y)); 
            graph_double.addVertex(ofVec2f(x,y+1)); 

            sum += U[n*nx + i];
        }
        graph.draw();
        graph_double.draw();

        ofSetColor(ofColor::red);
        ofDrawBitmapString("DTRW at time step " + ofToString(n+1) + " out of " + ofToString(nt) + " = " + ofToString(float(n) * dT, 2) + " seconds = " 
                           + ofToString(100.0 * float(n)/float(nt), 1) + "%", TW_MARGIN, TH_MARGIN);

        unlock();
    } 
}

void mcThread::start() 
{
    // Mutex blocking is set to true by default
    // It is rare that one would want to use startThread(false).
    startThread();
}
void mcThread::stop() 
{
        stopThread();    
}

void mcThread::threadedFunction() 
{
    while (isThreadRunning()) {
        if (lock()) {
            if (n < nMC-1) {
                n++;
                int t=0;
                int pos = nx / 2;
                t += inv_cum_dist(ran2(&seed), survival, nt);
                while (t < nt) {
                    if (ran2(&seed) > r) {
                        if (ran2(&seed) > 0.5) {
                            if (pos < nx-1) pos++;
                            //else pos--;
                        }
                        else { 
                            if (pos > 0) pos--;
                        //else pos++;
                        }
                    }
                    t += inv_cum_dist(ran2(&seed), survival, nt);
                }
                density[pos] += 1.0;
            }
            else
                stop();
            unlock();
        }

        else {
            // If we reach this else statement, it means that we could not
            // lock our mutex, and so we do not need to call unlock().
            // Calling unlock without locking will lead to problems.
            ofLogWarning("threadedFunction()") << "Unable to lock MC mutex.";
        }
    }
}
void mcThread::draw()
{
    if (lock()) {
        ofMesh graph;
        ofMesh graph_double;
        graph.setMode(OF_PRIMITIVE_LINE_STRIP);
        graph_double.setMode(OF_PRIMITIVE_LINE_STRIP);

        double sum = 0.0;
        for (int i=0; i<nx; i++) {
            float x = float((ofGetWidth()-2*W_MARGIN) * i) / float(nx) + float(W_MARGIN);
            float y = ofGetHeight() - ((ofGetHeight()-2*H_MARGIN) * (SCALE * density[i]/double(n)) + float(H_MARGIN));
            graph.addColor(ofColor::green);
            graph_double.addColor(ofColor::green);

            graph.addVertex(ofVec2f(x,y)); 
            graph_double.addVertex(ofVec2f(x,y+1)); 

            sum += double(density[i]) / double(n);
        }

        graph.draw();
        graph_double.draw();
    
        ofSetColor(ofColor::green);
        ofDrawBitmapString("MC at run " + ofToString(n+1) + " out of " + ofToString(nMC) + " = " + ofToString(100.0 * float(n) / float(nMC), 1) + "%", TW_MARGIN, TH_MARGIN + 12);

        unlock();
    }
}



void qmcThread::start() 
{
    // Mutex blocking is set to true by default
    // It is rare that one would want to use startThread(false).
    startThread();
}
void qmcThread::stop() 
{
        stopThread();    
}

void qmcThread::threadedFunction() 
{
    while (isThreadRunning()) {
        if (lock()) {
            if (n < nMC-1) {
                n++;

                int dim_count = 1;
                double* x_k = qmc_rule->next_point();
                double* y_k = Points::shift_points(x_k,shift,nt,false); // digital shift or normal shift
                int t=0;
                int pos = nx / 2;
                t += inv_cum_dist(y_k[dim_count], survival, nt);
                dim_count++;
                 
                while (t < nt) {
                    if (y_k[dim_count] > r) {
                        dim_count++;
                        if (y_k[dim_count] > 0.5) {
                            if (pos < nx-1) pos++;
                            //else pos--;
                        }
                        else { 
                            if (pos > 0) pos--;
                        //else pos++;
                        dim_count++;
                        }
                    }
                    t += inv_cum_dist(y_k[dim_count], survival, nt);
                    dim_count++;
                }
                density[pos] += 1.0;
                delete [] y_k;
                delete [] x_k;
            }
            else
                stop();
            unlock();
        }

        else {
            // If we reach this else statement, it means that we could not
            // lock our mutex, and so we do not need to call unlock().
            // Calling unlock without locking will lead to problems.
            ofLogWarning("threadedFunction()") << "Unable to lock MC mutex.";
        }
    }
}
void qmcThread::draw()
{
    if (lock()) {
        ofMesh graph;
        ofMesh graph_double;
        graph.setMode(OF_PRIMITIVE_LINE_STRIP);
        graph_double.setMode(OF_PRIMITIVE_LINE_STRIP);

        double sum = 0.0;
        for (int i=0; i<nx; i++) {
            float x = float((ofGetWidth()-2*W_MARGIN) * i) / float(nx) + float(W_MARGIN);
            float y = ofGetHeight() - ((ofGetHeight()-2*H_MARGIN) * (SCALE * density[i]/double(n)) + float(H_MARGIN));
            graph.addColor(ofColor::lightBlue);
            graph_double.addColor(ofColor::lightBlue);

            graph.addVertex(ofVec2f(x,y)); 
            graph_double.addVertex(ofVec2f(x,y+1)); 

            sum += double(density[i]) / double(n);
        }

        graph.draw();
        graph_double.draw();
    
        ofSetColor(ofColor::lightBlue);
        ofDrawBitmapString("QMC at run " + ofToString(n+1) + " out of " + ofToString(nMC) + " = " + ofToString(100.0 * float(n) / float(nMC), 1) + "%", TW_MARGIN, TH_MARGIN + 2*12);

        unlock();
    }
}

