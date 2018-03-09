#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
        return "";
    }
    else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

int main() {
    uWS::Hub h;

    // MPC is initialized here!
    MPC mpc;

    h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
        uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);

        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {
                    // j[1] is the data JSON object

                    ///////////////////////////////////
                    // Get Current State - where are we?
                    Dvector current_state(6);
                    Dvector current_actuators(2);
                    std::ofstream log;
                    
                    const double px = j[1]["x"];
                    const double py = j[1]["y"];
                    const double psi = j[1]["psi"];

                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];

                    if (LOG) {
                        log.open("logfile.txt",std::ios_base::app | std::ios_base::out);
                        log << "ptsx:" << "\n";
                        for (unsigned int i=0;i<ptsx.size();i++)
                            log << ptsx[i] << "\n";
                        log << "ptsy:" << "\n";
                        for (unsigned int i=0;i<ptsx.size();i++)
                            log << ptsy[i] << "\n";
                        log.close();
                    }
                    
                    for (unsigned int i = 0; i < ptsx.size(); i++) {

                        const double dx = ptsx[i] - px;
                        const double dy = ptsy[i] - py;

                        ptsx[i] = dx * cos(0-psi) - dy * sin(0-psi);
                        ptsy[i] = dy * cos(0-psi) + dx * sin(0-psi);
                    }

                    if (LOG) {
                        log.open("logfile.txt",std::ios_base::app | std::ios_base::out);
                        log << "px: " << px << "\n";
                        log << "py: " << py << "\n";
                        log << "psi: " << psi << "\n";
                        log << "ptsx transformed:" << "\n";
                        for (unsigned int i=0;i<ptsx.size();i++)
                            log << ptsx[i] << "\n";
                        log << "ptsy transformed:" << "\n";
                        for (unsigned int i=0;i<ptsx.size();i++)
                            log << ptsy[i] << "\n";
                        log.close();
                    }

                    ///////////////////////////////////
                    // Fit polynomial to waypoints ahead of us   
                    mpc.fitWaypoints(ptsx, ptsy);
                    if (DEBUG) {
                      for (auto i = ptsx.begin(); i != ptsx.end(); ++i)
                        std::cout << *i << ' ';
                    std::cout << std::endl;
                    for (auto i = ptsy.begin(); i != ptsy.end(); ++i)
                        std::cout << *i << ' ';                      
                    std::cout << std::endl;
                    std::cout << "coeffs" << mpc.coeffs << std::endl;
                }

                current_state[IDX_x] = 0;
                current_state[IDX_y] = 0;
                current_state[IDX_psi] = 0;                            
                    // Keep velocity v
                current_state[IDX_v] = j[1]["speed"];
                    // Estimate error at point (0,0)
                current_state[IDX_cte] = mpc.polyeval(current_state[IDX_x]) - current_state[IDX_y] ;
                    // epsi = psi - atan(coeffs[1] + 2 * px * coeffs[2] + 3 * coeffs[3] * pow(px, 2))
                current_state[IDX_epsi] = current_state[IDX_psi] - atan(mpc.polyeval_t(current_state[IDX_x]));                           

                current_actuators[IDX_delta] = j[1]["steering_angle"];
                current_actuators[IDX_a] = j[1]["throttle"];

                    ///////////////////////////////////
                    // GET THE CURRENT DELAYED STATE
                    // to cater for delayed response, update state to a state in 0.1s
                current_state = mpc.globalKinematic(current_state, current_actuators, 0.1);

                    // current state must be in vehicle coordinates with the delay factored in
                    // kinematic model is at play here
                    // note that at current state at vehicle coordinates:
                    // px, py, psi = 0.0, 0.0, 0.0
                    // note that in vehicle coordinates it is going straight ahead the x-axis
                    // which means position in vehicle's y-axis does not change
                    // the steering angle is negative the given value as we have
                    // as recall that during transformation we rotated all waypoints by -psi
                mpc.solve(current_state, current_actuators);

                json msgJson;
                msgJson["steering_angle"] = mpc.steer;
                msgJson["throttle"] = mpc.throttle;

                msgJson["next_x"] = mpc.next_xs;
                msgJson["next_y"] = mpc.next_ys;

                msgJson["mpc_x"] = mpc.mpc_xs;
                msgJson["mpc_y"] = mpc.mpc_ys;

                auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                if (DEBUG)
                    std::cout << msg << std::endl;

                this_thread::sleep_for(chrono::milliseconds(100));
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
        else {
                // Manual driving
            std::string msg = "42[\"manual\",{}]";
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
    }
});

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
    size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
        res->end(s.data(), s.length());
    }
    else {
            // i guess this should be done more gracefully?
        res->end(nullptr, 0);
    }
});

h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
});

h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
    char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
});

int port = 4567;
if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
}
else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
}
h.run();
}
