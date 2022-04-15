var layout = {
  title: "",
  xaxis: {
    tickangle: -90
  },
  yaxis: {
    title: {
      text: `${molecule} (${units})`,
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
