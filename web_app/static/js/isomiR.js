isomiR = window.location.href.split("/").pop();
isomiR = isomiR.split("#")[0];
d3.csv(`/api/isomiR/${isomiR}/expression`, function(err, rows){
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
    yaxis: {
      title: {
        text: `${isomiR} (RPM)`,
        standoff: 171717
      },
      zeroline: false
    },
    margin: {
      r: 0
    }
  }

  var config = {responsive: true}

  Plotly.newPlot("expression-plot", data, layout, config);
});
