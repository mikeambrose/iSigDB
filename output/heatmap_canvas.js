/*From http://www.highcharts.com/demo/heatmap-canvas */

$(function () {
xlabels = document.getElementById('xlabels').innerHTML.split(',');
ylabels = document.getElementById('ylabels').innerHTML.split(',');

//document.getElementById('container').css('cursor','crosshair');
    /**
     * This plugin extends Highcharts in two ways:
     * - Use HTML5 canvas instead of SVG for rendering of the heatmap squares. Canvas
     *   outperforms SVG when it comes to thousands of single shapes.
     * - Add a K-D-tree to find the nearest point on mouse move. Since we no longer have SVG shapes
     *   to capture mouseovers, we need another way of detecting hover points for the tooltip.
     */
    (function (H) {
        var Series = H.Series,
            each = H.each,
            wrap = H.wrap,
            seriesTypes = H.seriesTypes;

        /**
         * Create a hidden canvas to draw the graph on. The contents is later copied over 
         * to an SVG image element.
         */
        Series.prototype.getContext = function () {
            if (!this.canvas) {
                this.canvas = document.createElement('canvas');
                this.canvas.setAttribute('width', this.chart.chartWidth);
                this.canvas.setAttribute('height', this.chart.chartHeight);
                this.canvas.style.cursor = 'crosshair';
                this.canvas.style.position = 'absolute';
                this.image = this.chart.renderer.image('', 0, 0, this.chart.chartWidth, this.chart.chartHeight).add(this.group);

                this.ctx = this.canvas.getContext('2d');
            }
            return this.ctx;
        };

        /** 
         * Draw the canvas image inside an SVG image
         */ 
        Series.prototype.canvasToSVG = function () {
            this.image.attr({ href: this.canvas.toDataURL('image/png') });
        };

        /**
         * Wrap the drawPoints method to draw the points in canvas instead of the slower SVG,
         * that requires one shape each point.
         */
        H.wrap(H.seriesTypes.heatmap.prototype, 'drawPoints', function (proceed) {

            var ctx = this.getContext();
            
            if (ctx) {

                // draw the columns
                each(this.points, function (point) {
                    var plotY = point.plotY,
                        shapeArgs;

                    if (plotY !== undefined && !isNaN(plotY) && point.y !== null) {
                        shapeArgs = point.shapeArgs;

                        ctx.fillStyle = point.pointAttr[''].fill;
                        ctx.fillRect(shapeArgs.x, shapeArgs.y, shapeArgs.width, shapeArgs.height);
                    }
                });

                this.canvasToSVG();

            } else {
                this.chart.showLoading("Your browser doesn't support HTML5 canvas, <br>please use a modern browser");

                // Uncomment this to provide low-level (slow) support in oldIE. It will cause script errors on
                // charts with more than a few thousand points.
                //proceed.call(this);
            }
        });
    }(Highcharts));


    var start;
    $('#container').highcharts({

        data: {
            csv: document.getElementById('csv').innerHTML,
        },

        chart: {
            type: 'heatmap',
            height: 300+ylabels.length*20,
            width: 500+xlabels.length*20,
        },


        title: {
            text: document.getElementById('title').innerHTML,
            align: 'left',
            x: 40
        },

        subtitle: {
            text: null,
            align: 'left',
            x: 40
        },

        xAxis: {
            text: null,
            categories: xlabels,
            labels: {
               rotation: 90,
               style: {
                  fontSize: '10px',
                  step: 1
               }
            },
            //range from 0 to the number of elements
            tickPositions: Array.apply(null, Array(xlabels.length)).map(function (_, i) {return i;})
        },

        yAxis: {
            text: null,
            categories: ylabels,
        },
        legend: {
            title: {
                text: document.getElementById('legendLabel').innerHTML
            }
        },
            colorAxis: {
                min: parseFloat(document.getElementById('min').innerHTML),
                max: parseFloat(document.getElementById('max').innerHTML),
                stops: [
                    [0,'#00005C'],
                    [0.25, '#3060cf'],
                    [0.5, '#ffffff'],
                    [0.75, '#c4463a'],
                    [1,'#800000']
                ],
                //minColor: '#FFFFFF',
                //maxColor: Highcharts.getOptions().colors[0]
            },
       // },
	tooltip: {
            formatter: function () {
                return '<b>' + this.series.xAxis.categories[this.point.x] + ', ' +
                    this.series.yAxis.categories[this.point.y] + ', ' + this.point.value + '</b>';
            },
            //crosshairs: [true,true]
            },
        series: [{
            turboThreshold: 0
        }],

    });
//$('#window').resize();
    //$('#width').bind('input', function(){$('#container').highcharts().setSize(this.value,400,false);});
});
