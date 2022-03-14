//d3.csv("https://raw.githubusercontent.com/plotly/datasets/master/violin_data.csv", function(err, rows){
d3.csv("/api/isomiR/test/expression", function(err, rows){

  function unpack(rows, key) {
  return rows.map(function(row) { return row[key]; });
  }

var data = [{
  type: 'violin',
  y: unpack(rows, 'total_bill'),
  points: 'none',
  box: {
    visible: true
  },
  boxpoints: false,
  line: {
    color: 'black'
  },
  fillcolor: '#8dd3c7',
  opacity: 0.6,
  meanline: {
    visible: true
  },
  x0: "Total Bill"
}]

var layout = {
  title: "",
  yaxis: {
    zeroline: false
  }
}

Plotly.newPlot('myDiv', data, layout);
});
