molecule = window.location.href.split("/").pop();
molecule = molecule.split("#")[0];

var layout = {
  title: "",
  height: data.length * data[0]["num_cancers"] * 30,
  xaxis: {
    side: "top",
    title: {
      text: `${molecule} (RPM)`,
      standoff: 0
    },
    zeroline: false,
  },
  yaxis: {
  },
  margin: {
    r: 0
  },
  boxmode: "group"
}

var config = {responsive: true}

Plotly.newPlot("expression-plot", data, layout, config);
