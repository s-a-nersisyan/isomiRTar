molecule = window.location.href.split("/").pop();
molecule = molecule.split("#")[0];
d3.csv(`/api/molecule/${molecule}/expression`, function(err, rows){
  function unpack(rows, key) {
    return rows.map(function(row) { return row[key]; });
  }

  var data = [{
    type: "box",
    x: unpack(rows, "cancer"),
    y: unpack(rows, "expression")
  }]

  var layout = {
    title: "",
    xaxis: {
      tickangle: -90
    },
    yaxis: {
      title: {
        text: `${molecule} (RPM)`,
        standoff: 171717
      },
      zeroline: false,
    },
    margin: {
      r: 0
    }
  }

  var config = {responsive: true}

  Plotly.newPlot("expression-plot", data, layout, config);
});
