//
//  main.swift
//  baseball_path
//
//  Created by Wayne Cochran on 8/19/24.
//

import Foundation

let v0 = 208.0    // feet per sec
let angle = 43.0  // degrees
print("# v0 = \(v0), angle = \(angle)")

let path = BaseballPathGenerator.path2D(initialSpeedFeetPerSec: v0,
                                        initialAngleDegrees: angle)

print("# time = \(path.time)")
for point in path.path {
    print("\(point.x) \(point.y)")
}
