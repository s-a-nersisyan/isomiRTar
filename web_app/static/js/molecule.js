molecule = window.location.href.split("/").pop();
molecule = molecule.split("#")[0];
d3.json(`/api/molecule/${molecule}/expression`, function(err, data){

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