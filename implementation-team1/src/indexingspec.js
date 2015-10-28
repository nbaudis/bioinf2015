function alignedStair(n, align) {
  var repeats = Math.floor(n/align);
  var rest = n % align;
  var regularPart = repeats*(repeats+1)/2*align*align;
  var restPart = rest*(repeats + 1)*align;
  return regularPart + restPart;
}

function indexer(rows, cols, align) {
  return function(row, col) {
  
    var smallerRank = Math.min(rows, cols);
    var largerRank = Math.max(rows, cols);
    
    var diagonal = row + col;
    var d = diagonal;
    
    var closing = Math.max(d, largerRank) - largerRank; // essentially std::min(diagonal-largerRank, 0)
    d -= closing;

    // size of the middle part:
    var middle = Math.max(d, smallerRank) - smallerRank;
    d -= middle;

    // size of the opening part:
    var opening = d;

    // Now for computing the actual offsets. opening and middle part are straightforward:
    var offset = 0;
    
    offset += alignedStair(opening, align);
    console.log(offset);
    offset += diagonal < cols ? align - 1 : smallerRank % align == 1 ? 0 : align;

    offset += middle * Math.ceil(smallerRank/align)*align;

    // The closing part is more tricky, but upon rearranging I got to the following:
    offset += alignedStair(smallerRank - 1, align);
    offset -= alignedStair(smallerRank - 1 - closing, align);
    
    // What remains is indexing into the diagonal:
    offset += Math.min(row, cols - 1 - col);
    
    return offset;
  };
}


function check(rows, cols, align) {
  var idx = indexer(rows, cols, align);
  describe('indexer(' + rows + ', ' + cols + ', ' + align + ')', function() {
    var diagOffset = align - 1;
    for (var diag = 0; diag < rows + cols - 1; ++diag) {
      var top = Math.max(0, diag + 1 - cols);
      var bottom = Math.min(rows - 1, diag);
      var cells = bottom - top + 1;
      diagOffset += top == 0 
          ? (align - 1 - diagOffset % align) % align
          : (align - diagOffset % align) % align;
      
      describe('diagonal ' + diag, function() {
        for (var i = 0; i < 1; i++) {
          var row = top + i;
          var col = diag - top - i;
          var actual = idx(row, col);
          var expected = diagOffset + i;
          
          it('idx('+ row + ', ' + col + ')', function(){
            expect(actual).toBe(expected);
          });
        }
      });
      diagOffset += cells;
    }
  });
  
}

check(2, 4, 4);
check(7, 4, 4);
check(7, 5, 4);
check(7, 6, 4);
check(7, 7, 4);
check(7, 8, 4);
check(7, 9, 4);
check(7, 10, 4);
check(586, 556, 4);