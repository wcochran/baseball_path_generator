//
//  BaseballPathGenerator.swift
//  baseball_path
//
//  Created by Wayne Cochran on 8/19/24.
//

import Foundation
import simd

class BaseballPathGenerator {

    //
    // 4th order Runge-Kutta for system of N first-order differential equations.
    //  X' = F(t,X)
    //  X(a) = S
    // https://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
    //  S : initial value of X at t = t0
    //  t0 : initial time
    //  dt : time step
    //  maxIters : maximum number of steps
    //  terminate : function that determines (early) termination
    //  F : function that computes F(t,X).
    //
    private static func RK4(S : [Double], t0 : Double, dt : Double, maxIters : Int,
                            terminate : ([Double]) -> Bool,
                            F : (Double, [Double]) -> [Double]) -> (time: Double, samples: [[Double]]) {
        let N = S.count
        let n = maxIters
        let h = dt
        let h2 = h/2
        let h6 = h/6
        var X = S
        var samples : [[Double]] = [X]
        var tfinal = 0.0

        for k in 0 ..< n {
            let t = t0 + Double(k)*h
            let K1 = F(t, X)
            let Y1 = (0..<N).map {i -> Double in
                return X[i] + h2*K1[i]
            }
            let K2 = F(t+h2, Y1)
            let Y2 = (0..<N).map {i -> Double in
                return X[i] + h2*K2[i]
            }
            let K3 = F(t+h2, Y2)
            let Y3 = (0..<N).map {i -> Double in
                return X[i] + h*K3[i]
            }
            let K4 = F(t+h, Y3)
            for i in 0..<N {
                X[i] += h6*(K1[i] + 2*K2[i] + 2*K3[i] + K4[i])
            }

            if terminate(X) {
                break
            }

            samples.append(X)
            tfinal = t

        }

        return (time: tfinal + h, samples: samples)
    }

    private static func ballBelowGround(X : [Double]) -> Bool {
        assert(X.count == 4)
        return X[2] < 0
    }

    private static let g = 32.0  // acceleration due to gravity (ft/sec^2)

    //
    // Constant air resistance
    //
    private static func F_simple_air_resist(t : Double, X : [Double]) -> [Double] {
        assert(X.count == 4)
        let k = 0.25  // air resistance coefficient
        return [X[1], -k*X[1], X[3], -k*X[3] - g]
    }

    //
    // Air resistance proportional to square of the velocity of baseball
    //
    private static func F_air_resist_velocity_squared(t : Double, X : [Double]) -> [Double] {
        assert(X.count == 4)
        let k = 0.002
        let speed = sqrt(X[1]*X[1] + X[3]*X[3])
        return [X[1], -k*X[1]*speed, X[3], -k*X[3]*speed - g]
    }

    //
    // Compute baseball 2D path
    //
    static func path2D(initialSpeedFeetPerSec v0: Double,
                       initialAngleDegrees angle: Double,
                       initialHeightFeet y0 : Double = 3.0) -> (time: Double, path:[(x: Double, y: Double)]) {
        assert(v0 > 0)
        assert(y0 > 0)
        assert(0 < angle && angle < 90)

        //
        // Height of ball without wind resistance:
        // y(t) = -g/2*t^2 + v0*t + y0
        // Time of flight to reach ground at t where y(t) = 0
        //  t = (v0 + sqrt(v0^2 + 4*g/2*y0))/g
        //
        let vy0 = v0*sin(angle * Double.pi / 180)
        let extraTimeDueToAirResistance = 1.0 // upper bound
        let maxFlightTime = (vy0 + sqrt(vy0*vy0 + 2*g*y0))/g + extraTimeDueToAirResistance
        let dt = 1.0/64
        let n = Int(ceil(maxFlightTime/dt))
        let theta = angle * Double.pi / 180

        //
        // Solve for path use RK4 solver.
        //
        let S = [0, v0*cos(theta), y0, v0*sin(theta)]
        let soln = Self.RK4(S: S, t0: 0, dt: dt, maxIters: n,
                            terminate: Self.ballBelowGround,
                            // XXX F: F_simple_air_resist)
                            F : F_air_resist_velocity_squared)

        let path = soln.samples.map {X in
            return (x:X[0], y:X[2])
        }
        return (time: soln.time, path: path)
    }

    //
    // v0 : speed of ball off the bat (ft/sec)
    // verticalAngle : 0 < initial trajectory < 90
    // 0 <= horzAngle <= 90 for fairball (else foul ball)
    //
    static func path(initialSpeedFeetPerSec v0: Double,
                     initialVerticalAngleDegrees vertAngle: Double,
                     angleFromFirstBaseDegrees horzAngle: Double,
                     initialHeightFeet y0: Double = 3.0) -> (time: Double, path:[simd_float3]) {
        assert(v0 > 0)
        assert(y0 > 0)
        assert(0 < vertAngle && vertAngle < 90)

        let path2D = Self.path2D(initialSpeedFeetPerSec: v0,
                                 initialAngleDegrees: vertAngle,
                                 initialHeightFeet: y0)

        let theta = horzAngle * Double.pi / 180
        let cos_theta = cos(theta)
        let sin_theta = sin(theta)
        let path = path2D.path.map {coord -> simd_float3 in
            let x = coord.x * cos_theta
            let y = coord.x * sin_theta
            let z = coord.y
            return simd_float3(x: Float(x), y: Float(y), z: Float(z))
        }

        return (time: path2D.time, path: path)
    }

}
