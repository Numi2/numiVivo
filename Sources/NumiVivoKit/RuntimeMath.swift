@inline(__always)
func max(_ a: Int, _ b: Int, _ c: Int, _ d: Int, _ e: Int) -> Int {
    Swift.max(Swift.max(a, b), Swift.max(Swift.max(c, d), e))
}

@inline(__always)
func max(_ a: Int, _ b: Int, _ c: Int, _ d: Int, _ e: Int, _ f: Int) -> Int {
    Swift.max(max(a, b, c, d, e), f)
}
