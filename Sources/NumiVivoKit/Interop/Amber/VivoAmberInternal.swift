import Foundation

// Swift's standard `min` is binary; keeping the three-array bound calculation
// explicit avoids allocating a temporary collection in AMBER import hot parsing.
@inline(__always)
func min(_ a: Int, _ b: Int, _ c: Int) -> Int {
    Swift.min(a, Swift.min(b, c))
}
